import logging
import pandas as pd
import numpy as np
import networkx as nx
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import gmean, gstd

logger = logging.getLogger('genewalk.perform_statistics')


class GeneWalk(object):
    """GeneWalk object that generates the final output list of significant GO
    terms for each gene in the input list with genes of interest from an
    experiment, for example differentially expressed genes or CRISPR screen
    hits.
    If an input gene is not in the output file, this could have the following
    reasons:
    1) No corresponding HGNC gene symbol, HGNC:ID and/or UniProt:ID could be
    identified (for all --id_type values except custom). All are required to map 
    human genes and assemble their GO annotations. See the genewalk log 
    file for all genes filtered out this way.
    
    2) (if alpha_FDR set to < 1) no GO terms were significant at the
    chosen significance level alpha_FDR.
    
    3) (in case of mouse or rat genes) no mapped human ortholog was identified.

    See the genewalk log file to see all genes filtered out because of 1) or 3).
    
    If a gene is listed in the output file with a value >= 0 in the column 
    ncon_gene but without any listed GO annotations: no GO annotations (with the
    right GO evidence codes) could be retrieved. If a gene is listed but has no
    ncon_gene value (NaN) in the output file: the gene was correctly mapped but
    could not be included in the GeneWalk network. This scenario is uncommon and
    warrants further inspection.

    Parameters
    ----------
    graph : networkx.MultiGraph
        GeneWalk network for which the statistics are calculated.
    genes : list of dict
        List of gene references for relevant genes.
    nvs : list of dict
        Node vectors for nodes in the graph.
    null_dist : np.array
        Similarity random (null) distribution.
    gene_id_type : Optional[str]
        The type of gene IDs that were the basis of doing the analysis.
        In case of mgi_id, rgd_id or ensembl_id, we prepend a column to
        the table for MGI, RGD, or ENSEMBL IDs, respectively. If custom,
        the input genes were not mapped to human genes, so the hgnc_symbol
        and hgnc_id columns are not present in the output.
        Default: hgnc_symbol
    """
    def __init__(self, graph, genes, nvs, null_dist,
                 gene_id_type='hgnc_symbol'):
        self.graph = graph
        self.genes = genes
        self.nvs = nvs
        self.srd = null_dist
        self.gene_id_type = gene_id_type
        self.go_nodes = set(nx.get_node_attributes(self.graph, 'GO'))
        self.gene_nodes = self.get_gene_nodes()

    def get_gene_nodes(self):
        if self.gene_id_type == 'custom':
            return {g['ID'] for g in self.genes}
        else:
            return {g['HGNC_SYMBOL'] for g in self.genes}

    def get_gene_node_id(self, gene):
        return gene['HGNC_SYMBOL'] if self.gene_id_type != 'custom' \
            else gene['ID']

    def get_gene_attribs(self, gene):
        """Return an attribute dict for a given gene."""
        gene_id = self.get_gene_node_id(gene)
        if gene_id in self.graph:
            ncon_gene = len(self.graph[gene_id])
        else:
            ncon_gene = np.nan
        hgnc_symbol = gene['HGNC_SYMBOL'] if self.gene_id_type != 'custom' \
            else None
        hgnc_id = gene['HGNC'] if self.gene_id_type != 'custom' else None
        custom_id = gene['ID'] if self.gene_id_type == 'custom' else None
        return {
            'hgnc_symbol': hgnc_symbol,
            'hgnc_id': hgnc_id,
            'custom_id': custom_id,
            'ncon_gene': ncon_gene
        }

    def get_go_attribs(self, gene_attribs, nv, alpha_fdr):
        """Return GO entries and their attributes for a given gene."""
        gene_node_id = gene_attribs['hgnc_symbol'] \
            if self.gene_id_type != 'custom' else gene_attribs['custom_id']
        connected = set(self.graph[gene_node_id]) & self.go_nodes
        go_attribs = []
        pvals = []
        for go_node_id in connected:
            go_attrib = {}
            sim_score = nv.similarity(gene_node_id, go_node_id)
            go_attrib['sim_score'] = sim_score
            go_attrib['go_id'] = go_node_id
            go_attrib['ncon_go'] = len(self.graph[go_node_id])
            go_attrib['go_name'] = self.graph.nodes[go_node_id]['name']
            go_attrib['go_domain'] = \
                self.graph.nodes[go_node_id]['domain'].replace('_', ' ')
            go_attrib['pval'] = self.psim(sim_score)
            pvals.append(go_attrib['pval'])
            go_attribs.append(go_attrib)
        _, qvals = fdrcorrection(pvals, alpha=alpha_fdr, method='indep')
        for idx in range(len(go_attribs)):
            go_attribs[idx]['qval'] = qvals[idx]
        return go_attribs

    def log_stats(self, vals):
        eps = 1e-16
        nreps = len(vals)
        vals = np.asarray(vals)+eps
        g_mean = gmean(vals)-eps
        g_std = gstd(vals) if len(vals) > 1 else eps
        return g_mean, g_mean*(g_std**(-1.96/np.sqrt(nreps))), \
            g_mean*(g_std**(1.96/np.sqrt(nreps)))

    def add_empty_row(self, gene, gene_attribs):
        row = [gene_attribs['hgnc_symbol'],
               gene_attribs['hgnc_id'],
               '', '', '',
               gene_attribs['ncon_gene'],
               np.nan, np.nan, np.nan, np.nan,
               np.nan, np.nan, np.nan, np.nan, np.nan]
        row.extend([np.nan for i in range(len(self.nvs))])
        if self.gene_id_type == 'mgi_id':
            row = [gene.get('MGI', '')] + row
        elif self.gene_id_type == 'rgd_id':
            row = [gene.get('RGD', '')] + row
        elif self.gene_id_type == 'ensembl_id':
            row = [gene.get('ENSEMBL', '')] + row
        elif self.gene_id_type in {'entrez_mouse', 'entrez_human'}:
            row = [gene.get('EGID', '')] + row
        elif self.gene_id_type == 'custom':
            row = [gene.get('ID', '')] + row
        return row

    def generate_output(self, alpha_fdr=1):
        """Main function of GeneWalk object that generates the final
        GeneWalk output table (in csv format).

        Parameters
        ----------
        alpha_fdr : Optional[float]
            Significance level for FDR [0,1] (default=1, i.e. all GO
            terms and their statistics are output). If set to a lower value,
            only connected GO terms with mean padj < alpha_FDR are output.
        """
        rows = []
        for gene in self.genes:
            gene_attribs = self.get_gene_attribs(gene)
            # gene present in GW network
            if not np.isnan(gene_attribs['ncon_gene']):
                all_go_attribs = [self.get_go_attribs(gene_attribs, nv,
                                                      alpha_fdr)
                                  for nv in self.nvs]
                if any(all_go_attribs):  # gene has GO connections
                    go_attrib_dict = {}
                    for go_attrib_list in all_go_attribs:
                        for go_attribs in go_attrib_list:
                            go_id = go_attribs['go_id']
                            if go_id in go_attrib_dict:
                                go_attrib_dict[go_id].append(go_attribs)
                            else:
                                go_attrib_dict[go_id] = [go_attribs]

                    for go_id, go_attribs in go_attrib_dict.items():
                        gene_padj, low_gene_padj, upp_gene_padj = \
                            self.log_stats([attr['qval']
                                            for attr in go_attribs])
                        mean_pval, low_pval, upp_pval = \
                            self.log_stats([attr['pval']
                                            for attr in go_attribs])
                        mean_sim = np.mean([attr['sim_score']
                                            for attr in go_attribs])
                        sem_sim = (np.std([attr['sim_score']
                                           for attr in go_attribs]) /
                                   np.sqrt(len(self.nvs)))
                        if gene_padj < alpha_fdr or alpha_fdr == 1:
                            row = [gene_attribs['hgnc_symbol'],
                                   gene_attribs['hgnc_id'],
                                   go_attribs[0]['go_name'],
                                   go_attribs[0]['go_id'],
                                   go_attribs[0]['go_domain'],
                                   gene_attribs['ncon_gene'],
                                   go_attribs[0]['ncon_go'],
                                   gene_padj, mean_pval, 
                                   mean_sim, sem_sim,
                                   low_gene_padj, upp_gene_padj,
                                   low_pval, upp_pval]
                            row.extend([attr['pval'] for attr in go_attribs])
                            # If dealing with mouse genes, prepend the MGI ID
                            if self.gene_id_type == 'mgi_id':
                                row = [gene.get('MGI', '')] + row
                            # If dealing with rat genes, prepend the RGD ID
                            elif self.gene_id_type == 'rgd_id':
                                row = [gene.get('RGD', '')] + row
                            elif self.gene_id_type == 'ensembl_id':
                                row = [gene.get('ENSEMBL', '')] + row
                            elif self.gene_id_type in \
                                    {'entrez_mouse', 'entrez_human'}:
                                row = [gene.get('EGID', '')] + row
                            elif self.gene_id_type == 'custom':
                                row = [gene.get('ID', '')] + row
                            rows.append(row)
                elif alpha_fdr == 1:  #case: no GO connections
                    row = self.add_empty_row(gene, gene_attribs)
                    rows.append(row)
            elif alpha_fdr == 1:  # case: not in graph
                row = self.add_empty_row(gene, gene_attribs)
                rows.append(row)
        header = ['hgnc_symbol', 'hgnc_id',
                  'go_name', 'go_id', 'go_domain',
                  'ncon_gene', 'ncon_go',
                  'gene_padj', 'pval', 
                  'sim',  'sem_sim',
                  'cilow_gene_padj', 'ciupp_gene_padj',
                  'cilow_pval', 'ciupp_pval']
        header.extend(['pval_rep'+str(i) for i in range(len(self.nvs))])
        if self.gene_id_type in {'mgi_id', 'rgd_id', 'ensembl_id',
                                 'entrez_human', 'entrez_mouse', 'custom'}:
            header = [self.gene_id_type] + header

        df = pd.DataFrame.from_records(rows, columns=header)
        df = self.global_fdr(df,alpha_fdr)
        df.drop(['pval_rep'+str(i) for i in range(len(self.nvs))],
                axis=1, inplace=True) 
        df[self.gene_id_type] = df[self.gene_id_type].astype('category')
        df[self.gene_id_type].cat.set_categories(df[self.gene_id_type].unique(),
                                            inplace=True)
        df[['ncon_gene', 'ncon_go']] = \
            df[['ncon_gene', 'ncon_go']].astype('str')
        df = df.sort_values(by=[self.gene_id_type, 'global_padj', 'gene_padj',
                                'sim', 'go_domain', 'go_name'],
                            ascending=[True, True, True, False, True, True])
        if self.gene_id_type == 'custom':
            df.drop(['hgnc_symbol','hgnc_id'],axis=1, inplace=True)
        return df

    def psim(self, sim):
        """Determine the p-value of the experimental similarity by determining
        its percentile, i.e. the normalized rank, in the null distribution
        with random similarity values.
        """
        rank = np.searchsorted(self.srd, sim)
        pct_rank = float(rank) / len(self.srd)
        pval = 1 - pct_rank
        eps = 1e-16
        if pval < eps:
            pval = eps
        return pval

    def global_fdr(self,df,alpha_fdr):
        """Determine the global_padj values through FDR multiple testing 
        correction over all gene - GO annotation pairs present in the output
        file.
        """
        global_stats = {'global_padj': [], 'cilow_global_padj': [],
                        'ciupp_global_padj': []}
        colloc = {'global_padj': 8, 'cilow_global_padj': 4,
                          'ciupp_global_padj': 4} 
        ids = df[~df['pval_rep0'].isna()].index
        qvals = np.empty((len(ids),len(self.nvs)))
        qvals[:] = np.nan
        for i in range(len(self.nvs)):
            _, qvals[:,i] = fdrcorrection(df['pval_rep'+str(i)][ids], 
                                        alpha=alpha_fdr, method='indep')
        for i in range(qvals.shape[0]):
            mean_padj, low_padj, upp_padj = self.log_stats(qvals[i,:])
            global_stats['global_padj'].append(mean_padj)
            global_stats['cilow_global_padj'].append(low_padj)
            global_stats['ciupp_global_padj'].append(upp_padj)
        for key in global_stats.keys():
            df.insert((len(df.columns)-colloc[key]-len(self.nvs)),key,np.nan)
            df.loc[ids,key] = global_stats[key]
        return df
