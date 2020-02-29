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
    identified. All are required to map genes and assemble their GO
    annotations.  
    2) (if alpha_FDR set to < 1) no GO terms were significant at the
    chosen significance level alpha_FDR.  
    3) (in case of mouse genes) no mapped human ortholog was identified.  
    If a gene is listed in the output file with NaN values in the columns ncon_go 
    and ncon_gene, the gene was not included in the GeneWalk network. 

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
    """
    def __init__(self, graph, genes, nvs, null_dist):
        self.graph = graph
        self.genes = genes
        self.nvs = nvs
        self.srd = null_dist
        self.go_nodes = set(nx.get_node_attributes(self.graph, 'GO'))
        self.gene_nodes = set([g['HGNC_SYMBOL'] for g in self.genes])

    def get_gene_attribs(self, gene):
        """Return an attribute dict for a given gene."""
        if gene['HGNC_SYMBOL'] in self.graph:
            ncon_gene = len(self.graph[gene['HGNC_SYMBOL']])
        else:
            ncon_gene = np.nan
        return {
            'hgnc_symbol': gene['HGNC_SYMBOL'],
            'hgnc_id': gene['HGNC'],
            'ncon_gene': ncon_gene
        }

    def get_go_attribs(self, gene_attribs, nv, alpha_fdr):
        """Return GO entries and their attributes for a given gene."""
        gene_node_id = gene_attribs['hgnc_symbol']
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
        g_std = gstd(vals)
        return g_mean, g_mean*(g_std**(-1.96/np.sqrt(nreps))), \
            g_mean*(g_std**(1.96/np.sqrt(nreps)))

    def add_empty_row(self, gene, gene_attribs, base_id_type):
        row = [gene_attribs['hgnc_symbol'],
               gene_attribs['hgnc_id'],
               '', '', '',
               gene_attribs['ncon_gene'],
               np.nan, np.nan, np.nan, np.nan,
               np.nan, np.nan, np.nan, np.nan, np.nan]
        if base_id_type == 'mgi_id':
            row = [gene.get('MGI', '')] + row
        elif base_id_type == 'ensembl_id':
            row = [gene.get('ENSEMBL', '')] + row
        elif base_id_type == 'entrez_mouse' or \
             base_id_type == 'entrez_human':
            row = [gene.get('EGID', '')] + row
        return row

    def generate_output(self, alpha_fdr=1, base_id_type='hgnc_symbol'):
        """Main function of GeneWalk object that generates the final
        GeneWalk output table (in csv format).

        Parameters
        ----------
        alpha_fdr : Optional[float]
            Significance level for FDR [0,1] (default=1, i.e. all GO
            terms and their statistics are output). If set to a lower value,
            only connected GO terms with mean padj < alpha_FDR are output.
        base_id_type : Optional[str]
            The type of gene IDs that were the basis of doing the analysis.
            In case of mgi_id or ensembl_id, we prepend a column to the table
            for MGI or ENSEMBL IDs, respectively.
            Default: hgnc_symbol
        """
        rows = []
        for gene in self.genes:
            gene_attribs = self.get_gene_attribs(gene)
            # gene present in GW network
            if not np.isnan(gene_attribs['ncon_gene']):
                all_go_attribs = [self.get_go_attribs(gene_attribs, nv,
                                                      alpha_fdr)
                                  for nv in self.nvs]
                if all_go_attribs:  # gene has GO connections
                    go_attrib_dict = {}
                    for go_attrib_list in all_go_attribs:
                        for go_attribs in go_attrib_list:
                            go_id = go_attribs['go_id']
                            if go_id in go_attrib_dict:
                                go_attrib_dict[go_id].append(go_attribs)
                            else:
                                go_attrib_dict[go_id] = [go_attribs]

                    for go_id, go_attribs in go_attrib_dict.items():
                        mean_padj, low_padj, upp_padj = \
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
                        if mean_padj < alpha_fdr or alpha_fdr == 1:
                            row = [gene_attribs['hgnc_symbol'],
                                   gene_attribs['hgnc_id'],
                                   go_attribs[0]['go_name'],
                                   go_attribs[0]['go_id'],
                                   go_attribs[0]['go_domain'],
                                   gene_attribs['ncon_gene'],
                                   go_attribs[0]['ncon_go'],
                                   mean_padj, low_padj, upp_padj,
                                   mean_pval, low_pval, upp_pval,
                                   mean_sim, sem_sim]
                            # If dealing with mouse genes, prepend the MGI ID
                            if base_id_type == 'mgi_id':
                                row = [gene.get('MGI', '')] + row
                            elif base_id_type == 'ensembl_id':
                                row = [gene.get('ENSEMBL', '')] + row
                            elif base_id_type == 'entrez_mouse' or \
                                 base_id_type == 'entrez_human':
                                row = [gene.get('EGID', '')] + row
                            rows.append(row)
                elif alpha_fdr == 1:  # case: no GO connections
                    row = self.add_empty_row(gene, gene_attribs, base_id_type)
                    rows.append(row)
            elif alpha_fdr == 1:  # case: not in graph
                row = self.add_empty_row(gene, gene_attribs, base_id_type)
                rows.append(row)
        header = ['hgnc_symbol', 'hgnc_id',
                  'go_name', 'go_id', 'go_domain',
                  'ncon_gene', 'ncon_go',
                  'mean_padj', 'cilow_padj', 'ciupp_padj',
                  'mean_pval', 'cilow_pval', 'ciupp_pval',
                  'mean_sim',  'sem_sim',
                  ]
        if base_id_type in {'mgi_id', 'ensembl_id', 'entrez_human',
                            'entrez_mouse'}:
            header = [base_id_type] + header

        df = pd.DataFrame.from_records(rows, columns=header)
        df[base_id_type] = df[base_id_type].astype('category')
        df[base_id_type].cat.set_categories(df[base_id_type].unique(),
                                            inplace=True)
        df[['ncon_gene', 'ncon_go']] = \
            df[['ncon_gene', 'ncon_go']].astype('str')
        df = df.sort_values(by=[base_id_type, 'go_domain', 'mean_padj',
                                'mean_sim', 'go_name'],
                            ascending=[True, True, True, False, True])
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
