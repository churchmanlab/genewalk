import pandas as pd
import numpy as np
import pickle as pkl
import networkx as nx
from indra.databases import hgnc_client
import statsmodels.stats.multitest
from genewalk.genewalk.get_indra_stmts import load_genes
from genewalk.genewalk.nx_mg_assembler import Nx_MG_Assembler

class GeneWalk(object):
    """GeneWalk object that generates the final output list of significant GO terms
    for each gene in the input list with genes of interest from an experiment, eg DE genes or CRISPR screen hits. 
    If an input gene is not in the output, check GeneWalk_allGO.csv to see if it is present there:
    1) if not: there are no GO annotations or INDRA statements for this gene or 
    2) if present: no annotated GO terms were significant at the chosen significance level (alpha_FDR).

    Parameters
    ----------
    path : directory where files are generated (default '~/genewalk/'),
    path_GO : directory where GO ontology, of GOAtools are located (default '~/genewalk/GO/'),
    fgeneid : filename of input list with HGNC ids from genes of interest, (default: 'HGNCidForINDRA.csv'), 
    fstmts : pickle file with INDRA statements as generated with get_indra_stmts.py (default: 'HGNCidForINDRA.pkl'),
    fmg : pickle file with networkx multigraph as generated with nx_mg_assembler.py (default: 'GeneWalk_MG.pkl'),
    fnv : pickle file with node vectors as generated with deepwalk.py (default: 'GeneWalk_DW_nv.pkl'),
    fnull_dist : pickle file with null distributions for significance testing as generated
                with get_null_distributions.py (default: 'GeneWalk_DW_rand_simdists.pkl') 
    
    Attributes
    ----------
    hgncid : list of HGNC ids from genes of interest (loaded from fgeneid)
    MG : Nx_MG_Assembler object, with indra statements as MG.stmts (loaded from fstmts) and MG.graph (loaded from fmg)
    nv : node vectors (loaded from fnv) 
    srd : similarity random (null) distributions (loaded from fnull_dist)
    outdf : pandas.DataFrame that will be the final result of GeneWalk
    """
    
    def __init__(self,path='~/genewalk/',
                 fgeneid='HGNCidForINDRA.csv',
                 fstmts='HGNCidForINDRA.pkl',
                 fmg='GeneWalk_MG.pkl',
                 fnv='GeneWalk_DW_nv.pkl',
                 fnull_dist='GeneWalk_DW_rand_simdists.pkl',
                 path_GO='~/genewalk/GO/',
                 mouse_genes=False):
        self.path=path
        self.mouse_genes=mouse_genes
        if self.mouse_genes:
            self.mdf=pd.DataFrame()
            self._load_mouse_genes(fgeneid)#read mgi csv into self.mdf and mapped HGNC:ID
        else:
            self.hgncid=load_genes(self.path+fgeneid)#read hgnc list of interest
        self.outdf=pd.DataFrame()
        # Open pickled statements and initialize Nx_MG_Assembler
        with open(self.path+fstmts, 'rb') as f:
            stmts=pkl.load(f)
        self.MG=Nx_MG_Assembler(stmts,path_GO)
        del(stmts)    
        #load multigraph
        with open(self.path+fmg, 'rb') as f:
            self.MG.graph=pkl.load(f)
        self.GO_nodes=set(nx.get_node_attributes(self.MG.graph,'GO'))
        # load all node vectors    
        with open(self.path+fnv, 'rb') as f:
            self.nv=pkl.load(f)
        # Load similarity null distributions for significance testing       
        with open(self.path+fnull_dist, 'rb') as f:
            self.srd = pkl.load(f)
    
    def _load_mouse_genes(self,fname):
        """Return a list of human genes based on a table of mouse genes."""
        self.mdf = pd.read_csv(self.path+fname)
        mgi_ids = self.mdf['MGI Gene/Marker ID']
        genes = []
        for mgi_id in mgi_ids:
            if mgi_id.startswith('MGI:'):
                mgi_id = mgi_id[4:]
            hgnc_id = hgnc_client.get_hgnc_from_mouse(mgi_id)
            if not hgnc_id:
                hgnc_id = 'NA'
            genes.append(hgnc_id)
        self.mdf.insert(loc=0,column='HGNC:ID', value=pd.Series(genes, index=self.mdf.index))
    
    def generate_output(self,alpha_FDR=0.05,fname_out='GeneWalk.csv'): 
        """main function of GeneWalk object that generates the final output list 
        Parameters
        ----------
        alpha_FDR :  significance level for FDR (default=0.05)
        fname_out : filename of GeneWalk output file (default=GeneWalk.csv)
        """
        g_view=nx.nodes(self.MG.graph)
        if self.mouse_genes:
            self.outdf=pd.DataFrame(columns=['MGI','Symbol','mapped HGNC:ID','mapped HUGO','GO description','GO:ID',
                                                'N_con(gene)','N_con(GO)',
                                                'similarity','pval','padj'])
            hgncid=list(self.mdf['HGNC:ID'])
            for n in g_view:
                try: 
                    if self.MG.graph.node[n]['HGNC'] in hgncid:
                        hid=self.MG.graph.node[n]['HGNC']
                        mgis=self.mdf[self.mdf['HGNC:ID']==hid]['MGI Gene/Marker ID'].unique()
                        symbols=self.mdf[self.mdf['HGNC:ID']==hid]['Symbol'].unique()
                        N_gene_con=len(self.MG.graph[n])
                        for i in range(len(mgis)):
                            GOdf=self.get_GO_df(n,N_gene_con,alpha_FDR)
                            GOdf.insert(loc=0,column='MGI', 
                                        value=pd.Series(mgis[i], index=GOdf.index))
                            GOdf.insert(loc=1,column='Symbol', 
                                        value=pd.Series(symbols[i], index=GOdf.index))
                            GOdf.insert(loc=2,column='mapped HGNC:ID', 
                                        value=pd.Series(hid, index=GOdf.index))
                            GOdf.insert(loc=3,column='mapped HUGO', value=pd.Series(n, index=GOdf.index))
                            GOdf.insert(loc=6,column='N_con(gene)', value=pd.Series(N_gene_con, index=GOdf.index))
                            self.outdf=self.outdf.append(GOdf, ignore_index=True)
                except KeyError:
                    pass
            self.outdf=self.outdf.sort_values(by=['Symbol','padj','pval'])
        else:#human genes
            self.outdf=pd.DataFrame(columns=['HGNC:ID','HUGO','GO description','GO:ID',
                                                'N_con(gene)','N_con(GO)',
                                                'similarity','pval','padj'])
            for n in g_view:
                try: 
                    if self.MG.graph.node[n]['HGNC'] in self.hgncid:
                        N_gene_con=len(self.MG.graph[n])
                        GOdf=self.get_GO_df(n,N_gene_con,alpha_FDR)
                        GOdf.insert(loc=0,column='HGNC:ID', value=pd.Series(self.MG.graph.node[n]['HGNC'], index=GOdf.index))
                        GOdf.insert(loc=1,column='HUGO', value=pd.Series(n, index=GOdf.index))
                        GOdf.insert(loc=4,column='N_con(gene)', value=pd.Series(N_gene_con, index=GOdf.index))
                        self.outdf=self.outdf.append(GOdf, ignore_index=True)
                except KeyError:
                    pass
            self.outdf['HGNC:ID'] = self.outdf['HGNC:ID'].astype("category")
            self.outdf['HGNC:ID'].cat.set_categories(self.hgncid, inplace=True)#problematic if non-unique list self.hgncid
            self.outdf=self.outdf.sort_values(by=['HGNC:ID','padj','pval'])
        self.outdf.to_csv(self.path+fname_out, index=False)
        return self.outdf
    
    def P_sim(self,sim,N_con):
        dist_key='d'+str(np.floor(np.log2(N_con)))
        RANK = np.searchsorted(self.srd[dist_key], sim, side='left', sorter=None)
        PCT_RANK = float(RANK)/len(self.srd[dist_key])
        pval=1-PCT_RANK
        return pval

    def get_GO_df(self,geneoi,N_gene_con,alpha_FDR):
        N_GO_CON=[]
        PVAL=[]
        FDR=[]
        DES=[]
        GO_con2gene=set(self.MG.graph[geneoi]).intersection(self.GO_nodes)
        simdf=pd.DataFrame(self.nv.most_similar(geneoi,topn=len(self.nv.vocab)),columns=['GO:ID','similarity'])
        simdf=simdf[simdf['GO:ID'].isin(GO_con2gene)]
        
        for i in simdf.index:
            N_GO_con=len(self.MG.graph[simdf['GO:ID'][i]])
            N_GO_CON.append(N_GO_con)
            DES.append(self.MG.graph.node[simdf['GO:ID'][i]]['name'])
            pval=self.P_sim(simdf['similarity'][i],min(N_GO_con,N_gene_con))
            PVAL.append(pval)
        simdf.insert(loc=0,column='GO description', value=pd.Series(DES, index=simdf.index))
        simdf.insert(loc=2,column='N_con(GO)', value=pd.Series(N_GO_CON, index=simdf.index))
        simdf.insert(loc=4,column='pval', value=pd.Series(PVAL, index=simdf.index))
        BOOL,q_val=statsmodels.stats.multitest.fdrcorrection(simdf['pval'], 
                                                     alpha=alpha_FDR, method='indep')
        simdf.insert(loc=5,column='padj', value=pd.Series(q_val, index=simdf.index))
        return simdf[simdf['padj']<alpha_FDR]


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description='Choose a path where GeneWalk files are generated (default: ~/genewalk/ ).')
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--path_GO', default='(provided) path/GO/')
    parser.add_argument('--genes', default='data/JQ1_HGNCidForINDRA.csv')
    parser.add_argument('--stmts', default='data/JQ1_HGNCidForINDRA_stmts.pkl')
    parser.add_argument('--alpha_FDR', default=0.05)
    parser.add_argument('--mouse_genes',default=False)
    args = parser.parse_args()
    if args.path_GO=='(provided) path/GO/':
        args.path_GO=args.path+'GO/'

    GW=GeneWalk(path=args.path,
                    path_GO=args.path_GO,
                    fgeneid=args.genes,
                    fstmts=args.stmts,
                    fmg='GeneWalk_MG.pkl',
                    fnv='GeneWalk_DW_nv.pkl',
                    fnull_dist='GeneWalk_DW_rand_simdists.pkl',
                    mouse_genes=args.mouse_genes)
    GW.generate_output(alpha_FDR=args.alpha_FDR,fname_out='GeneWalk.csv')
    GW.generate_output(alpha_FDR=1,fname_out='GeneWalk_all_GO.csv')
