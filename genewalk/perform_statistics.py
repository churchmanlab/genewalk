import os
import copy
import logging
import argparse
import pandas as pd
import numpy as np
import pickle as pkl
import networkx as nx
import statsmodels.stats.multitest
from indra.databases import hgnc_client
from genewalk.get_indra_stmts import load_genes

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
    TODO: complete parameter docs

    Attributes
    ----------
    hgncid : list of str
        list of HGNC ids from genes of interest (loaded from fgenes in
        ase mouse_genes equals False)
    mdf : pandas.DataFrame
        pandas dataframe with MGI ids from genes of interest (loaded from
        fgenes in case mouse_genes equals True)
    graph : networkx.MultiGraph
        GeneWalk Network
    nv : dict
        node vectors (loaded from fnv_prefix and Nreps)
    srd : dict
        similarity random (null) distributions (loaded from fnull_dist)
    outdfs : list of pandas.DataFrame
        pandas DataFrames that will generate the final result of GeneWalk
    """

    # TODO: mouse gene mapping are loaded here to enable outputting the MGI
    #  IDs and symbols for mouse genes. This could be refactored to use INDRA's
    #  mappings and perhaps structured better in the code.
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
            go_attrib['ncon_go'] = self.graph[go_node_id]
            go_attrib['go_name'] = self.graph[go_node_id]['name']
            go_attrib['pval'] = \
                self.P_sim(sim_score, min(go_attrib['ncon_go'],
                                          gene_attribs['ncon_gene']))
            _, go_attrib['qval'] = \
                statsmodels.stats.multitest.fdrcorrection(go_attrib['pval'],
                                                          self.alpha_FDR,
                                                          method='indep')
            go_attribs.append(go_attrib)
        return go_attribs


    def generate_output(self, alpha_FDR=1):
        """Main function of GeneWalk object that generates the final output
        list

        Parameters
        ----------
        alpha_FDR
            significance level for FDR [0,1] (default=1, i.e. all GO
            terms are output). If set to a lower value, only annotated GO
            terms with mean padj < alpha_FDR are output.
        """
        for nv in self.nvs:
            g_view = nx.nodes(self.graph)
            for node in nx.nodes(self.graph):
                if self.graph.node[n]['HGNC'] in self.gene_nodes:
                    gene_attribs = self.get_gene_attribs(node)
                    go_attribs = self.get_go_attribs(gene_attribs, nv)

        #Merge all self.Nrep experimentals to calculate mean and sem statistics
        self.outdfs[self.Nreps+1]=copy.deepcopy(self.outdfs[1])
        for rep in range(2,self.Nreps+1):
            if self.mouse_genes: 
                COLUMNS=['MGI','Symbol','mapped HGNC','mapped Symbol',
                        'GO description','GO:ID',
                        'N_con(gene)','N_con(GO)']
            else:#human genes
                COLUMNS=['HGNC','Symbol',
                        'GO description','GO:ID',
                        'N_con(gene)','N_con(GO)']
            self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].merge(self.outdfs[rep],on=COLUMNS,how='outer')
        self.outdfs[self.Nreps+1].insert(loc=len(COLUMNS),column='mean:sim',
                    value=self.outdfs[self.Nreps+1][[str(r)+':similarity' for r in range(1,self.Nreps+1)]].mean(axis=1))
        self.outdfs[self.Nreps+1].insert(loc=len(COLUMNS)+1,column='sem:sim',
                    value=self.outdfs[self.Nreps+1][[str(r)+':similarity' for r in \
                                                          range(1,self.Nreps+1)]].std(axis=1)/np.sqrt(self.Nreps))
        self.outdfs[self.Nreps+1].insert(loc=len(COLUMNS)+2,column='mean:pval',
                    value=self.outdfs[self.Nreps+1][[str(r)+':pval' for r in range(1,self.Nreps+1)]].mean(axis=1))
        self.outdfs[self.Nreps+1].insert(loc=len(COLUMNS)+3,column='sem:pval',
                    value=self.outdfs[self.Nreps+1][[str(r)+':pval' for r in \
                                                          range(1,self.Nreps+1)]].std(axis=1)/np.sqrt(self.Nreps))
        self.outdfs[self.Nreps+1].insert(loc=len(COLUMNS)+4,column='mean:padj',
                    value=self.outdfs[self.Nreps+1][[str(r)+':padj' for r in range(1,self.Nreps+1)]].mean(axis=1))
        self.outdfs[self.Nreps+1].insert(loc=len(COLUMNS)+5,column='sem:padj',
                    value=self.outdfs[self.Nreps+1][[str(r)+':padj' for r in \
                                                          range(1,self.Nreps+1)]].std(axis=1)/np.sqrt(self.Nreps))   
        if self.mouse_genes:
            self.outdfs[self.Nreps+1]['MGI'] = self.outdfs[self.Nreps+1]['MGI'].astype("category")
            self.outdfs[self.Nreps+1]['MGI'].cat.set_categories(self.mdf['MGI'].unique(), inplace=True)
            self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].sort_values(by=['MGI','Symbol',
                                                                                'mean:padj','GO description'])
        else:#human genes
            self.outdfs[self.Nreps+1]['HGNC'] = self.outdfs[self.Nreps+1]['HGNC'].astype("category")
            self.outdfs[self.Nreps+1]['HGNC'].cat.set_categories(pd.Series(self.hgncid).unique(), inplace=True)
            self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].sort_values(by=['HGNC','mean:padj','GO description'])    
        self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].drop([str(r)+':similarity' for r in range(1,self.Nreps+1)],axis=1)
        self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].drop([str(r)+':pval' for r in range(1,self.Nreps+1)],axis=1)
        self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].drop([str(r)+':padj' for r in range(1,self.Nreps+1)],axis=1)
        return self.outdfs[self.Nreps+1]

    def P_sim(self, sim, N_con):
        # Gets the p-value by comparing the experimental similarity value
        # to the null distribution.
        # TODO: is searchsorted the slow step here?
        dist_key = 'd'+str(np.floor(np.log2(N_con)))
        RANK = np.searchsorted(self.srd[dist_key], sim, side='left',
                               sorter=None)
        PCT_RANK = float(RANK)/len(self.srd[dist_key])
        pval = 1-PCT_RANK
        return pval

    def get_GO_df(self, geneoi, N_gene_con, alpha_FDR):
        N_GO_CON = []
        PVAL = []
        FDR = []
        DES = []
        GO_con2gene = set(self.graph[geneoi]).intersection(self.GO_nodes)
        simdf = pd.DataFrame(self.nv.most_similar(geneoi,
                                                  topn=len(self.nv.vocab)),
                             columns=['GO:ID','similarity'])
        simdf=simdf[simdf['GO:ID'].isin(GO_con2gene)]

        for i in simdf.index:
            N_GO_con = len(self.graph[simdf['GO:ID'][i]])
            N_GO_CON.append(N_GO_con)
            DES.append(self.graph.node[simdf['GO:ID'][i]]['name'])
            pval = self.P_sim(simdf['similarity'][i],min(N_GO_con,N_gene_con))
            PVAL.append(pval)
        simdf.insert(loc=0,column='GO description',
                     value=pd.Series(DES, index=simdf.index))
        simdf.insert(loc=2,column='N_con(GO)',
                     value=pd.Series(N_GO_CON, index=simdf.index))
        simdf.insert(loc=4,column='pval',
                     value=pd.Series(PVAL, index=simdf.index))
        BOOL,q_val = \
            statsmodels.stats.multitest.fdrcorrection(simdf['pval'],
                                                      alpha=alpha_FDR,
                                                      method='indep')
        simdf.insert(loc=5, column='padj',
                     value=pd.Series(q_val, index=simdf.index))
        if alpha_FDR < 1:
            return simdf[simdf['padj'] < alpha_FDR]
        else:
            return simdf


if __name__ == '__main__':
    # Handle command line arguments
    # TODO: implement CLI with documentation here
    parser = argparse.ArgumentParser(
        description='Choose a path to the gene list.')
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--genes', default='gene_list.csv')
    parser.add_argument('--alpha_FDR', default=1)
    parser.add_argument('--mouse_genes',default=False)
    parser.add_argument('--filename_out',default='GeneWalk.csv')
    args = parser.parse_args()
    log_handler = logging.FileHandler(os.path.join(args.path, 'LogErr',
                                                   '%s.log' % logger.name))
    logger.addHandler(log_handler)
    GW = GeneWalk(path=args.path, fgenes=args.genes,
                  mouse_genes=args.mouse_genes)
    GW.generate_output(alpha_FDR=args.alpha_FDR, fname_out=args.filename_out)
