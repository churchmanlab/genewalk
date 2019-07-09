import logging
import pandas as pd
import numpy as np
import networkx as nx
import statsmodels.stats.multitest


logger = logging.getLogger('genewalk.perform_statistics')


class GeneWalk(object):
    """GeneWalk object that generates the final output list of significant GO
    terms for each gene in the input list with genes of interest from an
    experiment, for example differentially expressed genes or CRISPR screen
    hits.
    If an input gene is not in the output file, this could have the following
    reasons:
    1) No reaction statements could be retrieved from the data_source selected
    in get_node_vectors.py (Pathway Commons, indra or fromUser).
    2) No connected (annotated) GO terms are present in the GeneWalk Network.
    3) (In case of mouse genes) no mapped human ortholog was identified.
    4) (if alpha_FDR set to < 1) no annotated GO term were significant at the
    chosen significance level alpha_FDR.

    Parameters
    ----------
    graph : networkx.MultiGraph
        GeneWalk Network
    genes : list of dict
        List of gene references for relevant genes
    nvs : list of dict
        Node vectors
    null_dist : dict
        Similarity random (null) distributions
    """
    def __init__(self, graph, genes, nvs, null_dist):
        self.graph = graph
        self.genes = genes
        self.nvs = nvs
        self.srd = null_dist
        self.go_nodes = set(nx.get_node_attributes(self.graph, 'GO'))
        self.gene_nodes = set([g['HGNC_ID'] for g in self.genes])

    def get_gene_attribs(self, gene):
        return {
            'hgnc_symbol': gene['HGNC_SYMBOL'],
            'hgnc_id': gene['HGNC_ID'],
            'ncon_gene': len(self.graph[gene['HGNC_ID']])
        }

    def get_go_attribs(self, gene_attribs, nv):
        gene_node_id = gene_attribs['HGNC_ID']
        connected = set(self.graph[gene_node_id]) & self.GO_nodes
        similar = nv.most_similar(gene_node_id, topn=len(nv.vocab))
        connected_similar = [s for s in similar if s[0] in connected]
        go_attribs = []
        for go_node_id, sim_score in connected_similar:
            go_attrib = {}
            go_attrib['sim_score'] = sim_score
            go_attrib['go_id'] = go_node_id
            go_attrib['ncon_go'] = self.graph[go_node_id]
            go_attrib['go_name'] = self.graph[go_node_id]['name']
            go_attrib['pval'] = \
                self.psim(sim_score, min(go_attrib['ncon_go'],
                                         gene_attribs['ncon_gene']))
            _, go_attrib['qval'] = \
                statsmodels.stats.multitest.fdrcorrection(go_attrib['pval'],
                                                          self.alpha_FDR,
                                                          method='indep')
            go_attribs.append(go_attrib)
        return go_attribs

    def generate_output(self, alpha_fdr=1):
        """Main function of GeneWalk object that generates the final output
        list

        Parameters
        ----------
        alpha_FDR
            significance level for FDR [0,1] (default=1, i.e. all GO
            terms are output). If set to a lower value, only annotated GO
            terms with mean padj < alpha_FDR are output.
        """
        rows = []
        for node in self.gene_nodes:
            gene_attribs = self.get_gene_attribs(node)
            all_go_attribs = [self.get_go_attribs(gene_attribs, nv) for nv in
                              self.nvs]
            go_attrib_dict = {}
            for go_attrib_list in all_go_attribs:
                for go_attribs in go_attrib_list:
                    if go_attribs['go_id'] in go_attrib_dict:
                        go_attrib_dict[go_attribs['go_id']].append(go_attribs)
                    else:
                        go_attrib_dict[go_attribs['go_id']] = [go_attribs]

            for go_id, go_attribs in go_attrib_dict.items():
                mean_sim = np.mean([attr['sim_score'] for attr in go_attribs])
                ste_sim = (np.std([attr['sim_score'] for attr in go_attribs]) /
                           np.sqrt(len(self.nvs)))
                mean_pval = np.mean([attr['pval'] for attr in go_attribs])
                ste_pval = (np.std([attr['pval'] for attr in go_attribs]) /
                            np.sqrt(len(self.nvs)))
                mean_qval = np.mean([attr['qval'] for attr in go_attribs])
                ste_qval = (np.std([attr['qval'] for attr in go_attribs]) /
                            np.sqrt(len(self.nvs)))
                if mean_qval > alpha_fdr:
                    continue
                row = [gene_attribs['hgnc_symbol'],
                       gene_attribs['hgnc_id'],
                       go_attribs['go_name'],
                       go_attribs['go_id'],
                       gene_attribs['ncon_gene'],
                       go_attribs['ncon_go'],
                       mean_sim, ste_sim,
                       mean_pval, ste_pval,
                       mean_qval, ste_qval]
                rows.append(row)
        header = ['hgnc_symbol', 'hgnc_id', 'go_name', 'go_id', 'ncon_gene',
                  'ncon_go', 'mean_sim', 'ste_sim', 'mean_pval', 'ste_pval',
                  'mean_qval', 'ste_qval']
        df = pd.DataFrame.from_records(rows, columns=header)
        return df

    def psim(self, sim, ncon):
        # Gets the p-value by comparing the experimental similarity value
        # to the null distribution.
        # TODO: is searchsorted the slow step here?
        dist_key = 'd' + str(np.floor(np.log2(ncon)))
        rank = np.searchsorted(self.srd[dist_key], sim, side='left',
                               sorter=None)
        pct_rank = float(rank) / len(self.srd[dist_key])
        pval = 1 - pct_rank
        return pval
