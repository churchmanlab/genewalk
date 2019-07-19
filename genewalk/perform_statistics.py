import logging
import pandas as pd
import numpy as np
import networkx as nx
from statsmodels.stats.multitest import fdrcorrection

logger = logging.getLogger('genewalk.perform_statistics')


class GeneWalk(object):
    """GeneWalk object that generates the final output list of significant GO
    terms for each gene in the input list with genes of interest from an
    experiment, for example differentially expressed genes or CRISPR screen
    hits.
    If an input gene is not in the output file, this could have the following
    reasons:
    1) No corresponding HGNC gene symbol, HGNC:ID and/or UniProt:ID could be 
    identified. All are required to map genes and assemble their GO annotations.
    2) (in case of mouse genes) no mapped human ortholog was identified.
    
    Parameters
    ----------
    graph : networkx.MultiGraph
        GeneWalk network for which the statistics are calculated.
    genes : list of dict
        List of gene references for relevant genes.
    nvs : list of dict
        Node vectors for nodes in the graph.
    null_dist : dict
        Similarity random (null) distributions.
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
        pvals=[]
        for go_node_id in connected:
            go_attrib = {}
            sim_score = nv.similarity(gene_node_id,go_node_id)
            go_attrib['sim_score'] = sim_score
            go_attrib['go_id'] = go_node_id
            go_attrib['ncon_go'] = len(self.graph[go_node_id])
            go_attrib['go_name'] = self.graph.nodes[go_node_id]['name']
            go_attrib['pval'] = \
                self.psim(sim_score, min(go_attrib['ncon_go'],
                                         gene_attribs['ncon_gene']))            
            pvals.append(go_attrib['pval'])
            go_attribs.append(go_attrib)       
        _, qvals = fdrcorrection(pvals,alpha=alpha_fdr,method='indep')
        for idx in range(len(go_attribs)):
            go_attribs[idx]['qval'] = qvals[idx] #append the qval
            
        return go_attribs

    def generate_output(self, alpha_fdr=0.05, base_id_type='hgnc_symbol'):
        """Main function of GeneWalk object that generates the final 
        GeneWalk output table (in csv format).

        Parameters
        ----------
        alpha_fdr : Optional[float]
            Significance level for FDR [0,1] (default = 0.05).
        base_id_type : Optional[str]
            The type of gene IDs that were the basis of doing the analysis.
            In case of mgi_id, we prepend a column to the table for MGI IDs.
            Default: hgnc_symbol
        """
        rows = []
        for gene in self.genes:
            gene_attribs = self.get_gene_attribs(gene)
            if gene_attribs['ncon_gene'] > 0:
                all_go_attribs = [self.get_go_attribs(gene_attribs, nv, alpha_fdr)
                                  for nv in self.nvs]
                if all_go_attribs:#gene has GO connections
                    go_attrib_dict = {}
                    for go_attrib_list in all_go_attribs:
                        for go_attribs in go_attrib_list:
                            if go_attribs['go_id'] in go_attrib_dict:
                                go_attrib_dict[go_attribs['go_id']].append(go_attribs)
                            else:
                                go_attrib_dict[go_attribs['go_id']] = [go_attribs]

                    for go_id, go_attribs in go_attrib_dict.items():
                        mean_padj = np.mean([attr['qval'] for attr in go_attribs])
                        sem_padj = (np.std([attr['qval'] for attr in go_attribs]) /
                                    np.sqrt(len(self.nvs)))
                        mean_sim = np.mean([attr['sim_score'] for attr in go_attribs])
                        sem_sim = (np.std([attr['sim_score'] for attr in go_attribs]) /
                                   np.sqrt(len(self.nvs)))
                        mean_pval = np.mean([attr['pval'] for attr in go_attribs])
                        sem_pval = (np.std([attr['pval'] for attr in go_attribs]) /
                                    np.sqrt(len(self.nvs)))
                        row = [gene_attribs['hgnc_symbol'],
                               gene_attribs['hgnc_id'],
                               go_attribs[0]['go_name'],
                               go_attribs[0]['go_id'],
                               gene_attribs['ncon_gene'],
                               go_attribs[0]['ncon_go'],
                               mean_padj, sem_padj,
                               mean_pval, sem_pval,
                               mean_sim, sem_sim,
                           ]
                        # If we're dealing with mouse genes, prepend the MGI ID
                        if base_id_type == 'mgi_id':
                            row = [gene.get('MGI', '')] + row
                        rows.append(row)
                else:#case: no GO connections
                    row = [gene_attribs['hgnc_symbol'],
                           gene_attribs['hgnc_id'],
                           '',
                           '',
                           gene_attribs['ncon_gene'],
                           len(all_go_attribs),
                           np.nan, np.nan,
                           np.nan, np.nan,
                           np.nan, np.nan,
                       ]
                    if base_id_type == 'mgi_id':
                        row = [gene.get('MGI', '')] + row
                    rows.append(row)
            else:#case: no connections or not in graph
                row = [gene_attribs['hgnc_symbol'],
                       gene_attribs['hgnc_id'],
                       '',
                       '',
                       gene_attribs['ncon_gene'],
                       np.nan,
                       np.nan, np.nan,
                       np.nan, np.nan,
                       np.nan, np.nan,
                    ]
                if base_id_type == 'mgi_id':
                    row = [gene.get('MGI', '')] + row
                rows.append(row)
        header = ['hgnc_symbol', 'hgnc_id',
                  'go_name', 'go_id',
                  'ncon_gene', 'ncon_go',
                  'mean_padj', 'sem_padj',
                  'mean_pval', 'sem_pval',
                  'mean_sim', 'sem_sim',
                  ]
        if base_id_type == 'mgi_id':
            header = ['mgi_id'] + header
        df = pd.DataFrame.from_records(rows, columns=header)
        df[base_id_type] = df[base_id_type].astype('category')
        df[base_id_type].cat.set_categories(df[base_id_type].unique(), inplace=True)
        #TODO: decide to require pandas v0.24 with installation of GeneWalk 
        df[['ncon_gene', 'ncon_go']] = df[['ncon_gene', 'ncon_go']].astype('str')
        #df[['ncon_gene', 'ncon_go']] = df[['ncon_gene', 'ncon_go']].astype(pd.Int32Dtype())
        #better but works only for pandas >= v0.24:important to output as integer (with nan values)
        df = df.sort_values(by=[base_id_type,'mean_padj','go_name']) 
        return df

    def psim(self, sim, ncon):
        """
        Determine the p-value of the experimental similarity by determining 
        its percentile, i.e. the normalized rank, in the null distribution with 
        random similarity values.
        """
        dist_key = 'd' + str(np.floor(np.log2(ncon)))
        rank = np.searchsorted(self.srd[dist_key], sim)
        pct_rank = float(rank) / len(self.srd[dist_key])
        pval = 1 - pct_rank
        return pval
