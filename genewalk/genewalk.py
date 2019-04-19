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
    1) if not: there are no GO annotations or INDRA statements for this gene (or in case of mouse: no mapped human ortholog) or 
    2) if present: no annotated GO terms were significant at the chosen significance level (alpha_FDR).

    Parameters
    ----------
    path : directory where files are generated (default '~/genewalk/'),
    path_GO : directory where GO ontology, of GOAtools are located (default '~/genewalk/GO/'),
    fgeneid : filename of input list with HGNC ids (or MGI ids) from genes of interest, (default: 'HGNCidForINDRA.csv'), 
    fstmts : pickle file with INDRA statements as generated with get_indra_stmts.py (default: 'HGNCidForINDRA.pkl'),
    fmg : pickle file with networkx multigraph as generated with nx_mg_assembler.py (default: 'GeneWalk_MG.pkl'),
    fnv_prefix : prefix for multiple pickle files with node vectors as generated with deepwalk.py (default: 'GeneWalk_DW_nv'),
    Nreps : number of DeepWalk repetitions for the gene/GO network. It must match Nreps parameter used in deepwalk.py
            (default: 10), 
    fnull_dist : pickle file with null distributions for significance testing as generated
                with get_null_distributions.py (default: 'GeneWalk_DW_rand_simdists.pkl')
    mouse_genes : set to True if the input list are MGI:IDs from mouse genes (default: False)
    
    Attributes
    ----------
    hgncid : list of HGNC ids from genes of interest (loaded from fgeneid in case mouse_genes equals False)
    mdf : pandas dataframe with MGI ids from genes of interest (loaded from fgeneid in case mouse_genes equals True)
    MG : Nx_MG_Assembler object, with indra statements as MG.stmts (loaded from fstmts) and MG.graph (loaded from fmg)
    nv : node vectors (loaded from fnv_prefix and Nreps) 
    srd : similarity random (null) distributions (loaded from fnull_dist)
    outdfs : pandas.DataFrames that will generate the final result of GeneWalk
    """
    
    def __init__(self,path='~/genewalk/',
                 fgeneid='HGNCidForINDRA.csv',
                 fstmts='HGNCidForINDRA.pkl',
                 fmg='GeneWalk_MG.pkl',
                 fnv_prefix='GeneWalk_DW_nv',
                 Nreps=10,
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
        # Open pickled statements and initialize Nx_MG_Assembler
        with open(self.path+fstmts, 'rb') as f:
            stmts=pkl.load(f)
        self.MG=Nx_MG_Assembler(stmts,path_GO)
        del(stmts)    
        #load multigraph
        with open(self.path+fmg, 'rb') as f:
            self.MG.graph=pkl.load(f)
        self.GO_nodes=set(nx.get_node_attributes(self.MG.graph,'GO'))
        self.fnv_prefix=fnv_prefix
        self.Nreps=Nreps
        self.nv=[]#node vectors, defined in generate_output
        # Load similarity null distributions for significance testing       
        with open(self.path+fnull_dist, 'rb') as f:
            self.srd = pkl.load(f)
        self.outdfs=dict()#pd.DataFrame()
    
    def _load_mouse_genes(self,fname):
        """Append human gene IDs to a df of mouse genes (self.mdf)."""
        self.mdf = pd.read_csv(self.path+fname)#assumes the csv has headers
        for c in self.mdf.columns:
            if c.startswith('MGI'):#assumes the first column starting with MGI is the relevant one with MGI:IDs
                self.mdf=self.mdf.rename(columns={c: 'MGI'})
                break
        mgi_ids = self.mdf['MGI']
        genes = []
        for mgi_id in mgi_ids:
            if mgi_id.startswith('MGI:'):
                mgi_id = mgi_id[4:]
            hgnc_id = hgnc_client.get_hgnc_from_mouse(mgi_id)
            if not hgnc_id:
                hgnc_id = 'NA'
            genes.append(hgnc_id)
        self.mdf.insert(loc=0,column='HGNC:ID', value=pd.Series(genes, index=self.mdf.index))
    
    def generate_output(self,alpha_FDR=1.001,fname_out='GeneWalk.csv'): 
        """main function of GeneWalk object that generates the final output list 
        Parameters
        ----------
        alpha_FDR : significance level for FDR [0,1] (default=1, i.e. all GO terms are output). 
                    If set to a lower value, only annotated GO terms with mean padj < alpha_FDR are output.
        fname_out : filename of GeneWalk output file (default=GeneWalk.csv)
        """
        if self.mouse_genes: 
            hgncid=list(self.mdf['HGNC:ID'])
        else:#human genes
            hgncid=self.hgncid    
        for rep in range(1,self.Nreps+1):
            print(rep,'/',self.Nreps)

            # load node vectors
            fnv=self.fnv_prefix+'_'+str(rep)+'.pkl'
            with open(self.path+fnv, 'rb') as f:
                self.nv=pkl.load(f)
            g_view=nx.nodes(self.MG.graph)
            
            # initialize GeneWalk output dataframe for each replicate run  
            if self.mouse_genes: 
                COLUMNS=['MGI','Symbol','mapped HGNC:ID','mapped HUGO',
                                                       'GO description','GO:ID',
                                                       'N_con(gene)','N_con(GO)',
                                                       'similarity','pval','padj']
            else:#human genes
                COLUMNS=['HGNC:ID','HUGO',
                           'GO description','GO:ID',
                           'N_con(gene)','N_con(GO)',
                           'similarity','pval','padj']
            self.outdfs[rep]=pd.DataFrame(columns=COLUMNS)    
            
            if self.mouse_genes:
                for n in g_view:
                    try: 
                        if self.MG.graph.node[n]['HGNC'] in hgncid:
                            hid=self.MG.graph.node[n]['HGNC']
                            mgis=self.mdf[self.mdf['HGNC:ID']==hid]['MGI'].unique()
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
                                self.outdfs[rep]=self.outdfs[rep].append(GOdf, ignore_index=True)
                    except KeyError:
                        pass
                self.outdfs[rep]['MGI'] = self.outdfs[rep]['MGI'].astype("category")
                self.outdfs[rep]['MGI'].cat.set_categories(self.mdf['MGI'].unique(), inplace=True)
                self.outdfs[rep]=self.outdfs[rep].sort_values(by=['MGI','mapped HGNC:ID','mapped HUGO','GO:ID'])
            else:#human genes
                for n in g_view:
                    try: 
                        if self.MG.graph.node[n]['HGNC'] in hgncid:
                            N_gene_con=len(self.MG.graph[n])
                            GOdf=self.get_GO_df(n,N_gene_con,alpha_FDR)
                            GOdf.insert(loc=0,column='HGNC:ID', value=pd.Series(self.MG.graph.node[n]['HGNC'],
                                                                                index=GOdf.index))
                            GOdf.insert(loc=1,column='HUGO', value=pd.Series(n, index=GOdf.index))
                            GOdf.insert(loc=4,column='N_con(gene)', value=pd.Series(N_gene_con,index=GOdf.index))
                            self.outdfs[rep]=self.outdfs[rep].append(GOdf, ignore_index=True)
                    except KeyError:
                        pass
                self.outdfs[rep]['HGNC:ID'] = self.outdfs[rep]['HGNC:ID'].astype("category")
                self.outdfs[rep]['HGNC:ID'].cat.set_categories(pd.Series(self.hgncid).unique(), inplace=True)
                self.outdfs[rep]=self.outdfs[rep].sort_values(by=['HGNC:ID','HUGO','GO:ID'])
        
            mppr={'similarity':str(rep)+':similarity','pval':str(rep)+':pval','padj':str(rep)+':padj'}
            self.outdfs[rep]=self.outdfs[rep].rename(mapper=mppr,axis=1)
        #Merge all self.Nrep experimentals to calculate mean and sem statistics
        self.outdfs[self.Nreps+1]=copy.deepcopy(self.outdfs[1])        
        for rep in range(2,self.Nreps+1):
            if self.mouse_genes: 
                COLUMNS=['MGI','Symbol','mapped HGNC:ID','mapped HUGO',
                        'GO description','GO:ID',
                        'N_con(gene)','N_con(GO)']
            else:#human genes
                COLUMNS=['HGNC:ID','HUGO',
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
            self.outdfs[self.Nreps+1]['HGNC:ID'] = self.outdfs[self.Nreps+1]['HGNC:ID'].astype("category")
            self.outdfs[self.Nreps+1]['HGNC:ID'].cat.set_categories(pd.Series(self.hgncid).unique(), inplace=True)
            self.outdfs[self.Nreps+1]=self.outdfs[self.Nreps+1].sort_values(by=['HGNC:ID','mean:padj','GO description'])    
        self.outdfs[self.Nreps+1].drop([str(r)+':similarity' for r in range(1,self.Nreps+1)],axis=1)
        self.outdfs[self.Nreps+1].drop([str(r)+':pval' for r in range(1,self.Nreps+1)],axis=1)
        self.outdfs[self.Nreps+1].drop([str(r)+':padj' for r in range(1,self.Nreps+1)],axis=1)    
        self.outdf[self.Nreps+1].to_csv(self.path+fname_out, index=False)
        return self.outdf[self.Nreps+1]
    
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
        description='Choose a path where GeneWalk files are generated (default: ~/genewalk/ ). \
        Nreps has to match what is used  in get_node_vectors.py')
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--path_GO', default='(provided) path/GO/')
    parser.add_argument('--genes', default='data/JQ1_HGNCidForINDRA.csv')
    parser.add_argument('--stmts', default='data/JQ1_HGNCidForINDRA_stmts.pkl')
    parser.add_argument('--alpha_FDR', default=1.001)
    parser.add_argument('--mouse_genes',default=False)
    parser.add_argument('--filename_out',default='GeneWalk.csv')
    parser.add_argument('--Nreps', default=10)
    args = parser.parse_args()
    if args.path_GO=='(provided) path/GO/':
        args.path_GO=args.path+'GO/'

    GW=GeneWalk(path=args.path,
                    path_GO=args.path_GO,
                    fgeneid=args.genes,
                    fstmts=args.stmts,
                    fmg='GeneWalk_MG.pkl',
                    fnv_prefix='GeneWalk_DW_nv',
                    Nrep=args.Nreps,
                    fnull_dist='GeneWalk_DW_rand_simdists.pkl',
                    mouse_genes=args.mouse_genes)
    GW.generate_output(alpha_FDR=args.alpha_FDR,fname_out=args.filename_out)
