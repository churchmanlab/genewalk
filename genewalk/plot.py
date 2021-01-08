import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

logger = logging.getLogger('genewalk.plot')

class GW_Plotter(object):
    """ 
    genewalk plotter object to visualize genewalk output data.
    
    Parameters
    ----------
    path : path to output figures (default: project_folder/figures)
    dgw : genewalk_results pandas.dataFrame
    alpha_fdr : FDR significance level for global analyses and barplots
                (default: 0.1)
    """
    def __init__(self, path, dgw, alpha_fdr):
        self.path = path
        if alpha_fdr == 1:#default input
            self.alpha_fdr = 0.1
        else:
            self.alpha_fdr = alpha_fdr
        self.scatter_data = []
        self.stat = 'gene_padj' #for barplots 
        self.ci_stat = 'cilow_'+self.stat
        self.dGW = dgw
        self.go_domains = set(self.dGW['go_domain'])
        if self.dGW.columns[0] in {'mgi_id', 'ensembl_id', 
                                   'entrez_human','entrez_mouse'}:
            self.std_id = False
            self.id_type = self.dGW.columns[0]
        else:
            self.std_id = True
            self.id_type = self.dGW.columns[1]
    
    
    def generate_plots(self):
        """
        Wrapper that calls scatter and bar plot
        generating functions.
        """
        self.scatterplot_regulators()
        self.scatterplot_moonlighters()
        self.barplot_goanno()
    
    
    def scatterplot_regulators(self):
        """Scatter plot with fraction of (globally) relevant GO annotations 
        as a function of gene connectivity (to other genes) for all input
        genes. Genes with symbols listed are regulator genes.
        See genewalk_scatterplots.csv for full data and 
        genewalk_regulators.csv for regulator genes of interest.
        Visualization thresholds:
        T_frac: minimal fraction of relevant GO annotations, set to 0.5
        T_gcon: minimal gene connectivity (to other 
        genes), set to 75th quartile of distribution
        """
        if not self.scatter_data:
            self._get_scatter_data()
        xvar = 'gene_con'
        yvar = 'frac_rel_go'
        T_gcon = self.scatter_data[xvar].quantile(q=0.75)
        T_frac = 0.5
        sns.set(style="whitegrid")    
        fig, ax = plt.subplots(figsize=(10,10))#inches
        g = sns.scatterplot(x=xvar, y=yvar, hue=yvar, 
                            linewidth=0, alpha=0.5,
                            sizes=(40, 400), 
                            data=self.scatter_data,
                            ax=ax,legend=False)
        plt.axvline(x=T_gcon, color=[0.7,0.7,0.7], linestyle='--')
        plt.axhline(y=T_frac, color=[0.7,0.7,0.7], linestyle='--')
        font_sz=16
        plt.xlabel(xvar,size=font_sz)
        plt.ylabel(yvar,size=font_sz)
        plt.xlim([0,max(self.scatter_data[xvar])*1.2])
        plt.xticks(size=font_sz)
        plt.yticks(size=font_sz)

        regulators=[]
        dreg = self.scatter_data[ self.scatter_data[xvar] >= T_gcon ]
        dreg = dreg[ dreg[yvar] >= T_frac ]
        for r in dreg.index:
            gname = dreg['hgnc_symbol'][r]
            regulator.append(gname)  
            x_txt=dreg[xvar][t]
            y_txt=dreg[yvar][t]#+0.01
            g.text(x_txt, y_txt, gname, size=6,horizontalalignment='center', 
                   color='black',weight='semibold',
                   fontstyle='italic')
        g.set(xscale="log")
        plt.title('Regulator genes',size=font_sz)
        filename = 'regulators_x_'+xvar+'_y_'+yvar+'.'
        plt.savefig(os.path.join(self.path,filename+'pdf'),
                    bbox_inches="tight",transparent=True)
        plt.savefig(os.path.join(self.path,filename+'png'),
                    bbox_inches="tight",transparent=True)
        
        df = pd.DataFrame(regulators, columns=['gw_regulators'])
        df.to_csv(os.path.join(self.path,'genewalk_regulators.csv'),index=False)
    
    
    def scatterplot_moonlighters(self):
        """Scatter plot with fraction of (globally) relevant GO annotations 
        as a function of number of GO annotations for all input genes.
        Genes with symbols listed are moonlighting genes.
        See genewalk_scatterplots.csv for full data and 
        genewalk_moonlighters.csv for moonlighting genes of interest.
        Visualization thresholds:
        T_frac: maximal fraction of relevant GO annotations, set to 0.5
        T_gocon: minimal number of GO annotations per gene, 
        set to max of 30 and 75th quartile of distribution.
        """
        if not self.scatter_data:
            self._get_scatter_data()    
        xvar = 'go_con'
        yvar = 'frac_rel_go'
        T_gocon = max(30,self.scatter_data[xvar].quantile(q=0.75))
        T_frac = 0.5
        sns.set(style="whitegrid")    
        fig, ax = plt.subplots(figsize=(10,10))#inches
        g = sns.scatterplot(x=xvar, y=yvar, hue=yvar, 
                            linewidth=0, alpha=0.5,
                            sizes=(40, 400), 
                            data=self.scatter_data,
                            ax=ax,legend=False)
        plt.axvline(x=T_gocon, color=[0.7,0.7,0.7], linestyle='--')
        plt.axhline(y=T_frac, color=[0.7,0.7,0.7], linestyle='--')
        font_sz=16
        plt.xlabel(xvar,size=font_sz)
        plt.ylabel(yvar,size=font_sz)
        plt.xlim([0,max(self.scatter_data[xvar])*1.2])
        plt.xticks(size=font_sz)
        plt.yticks(size=font_sz)

        moonlighters=[]
        dmoon = self.scatter_data[ self.scatter_data[xvar] >= T_gocon ]
        dmoon = dmoon[ dmoon[yvar] <= T_frac ]
        for r in dreg.index:
            gname = dreg['hgnc_symbol'][r]
            moonlighters.append(gname)  
            x_txt=dreg[xvar][t]
            y_txt=dreg[yvar][t]#+0.01
            g.text(x_txt, y_txt, gname, size=6,horizontalalignment='center', 
                   color='black',weight='semibold',
                   fontstyle='italic')
        g.set(xscale="log")
        plt.title('Moonlighting genes',size=font_sz)
        filename = 'moonlighters_x_'+xvar+'_y_'+yvar+'.'
        plt.savefig(os.path.join(self.path,filename+'pdf'),
                    bbox_inches="tight",transparent=True)
        plt.savefig(os.path.join(self.path,filename+'png'),
                    bbox_inches="tight",transparent=True)
        
        df = pd.DataFrame(moonlighters, columns=['gw_moonlighter'])
        df.to_csv(os.path.join(self.path,'genewalk_moonlighters.csv'),index=False)
    
    
    def barplot_goanno(self):
        """Visualize statistical significances of GO annotations for a given 
        gene of interest. Four separate plots are generated: one with all GO 
        annotations listed and 3 separated by each go domain: biological process, 
        cellular component and molecular function.
        """ 
        self.dGW['mlog10padj'] = -np.log10(self.dGW['gene_padj'])
        self.dGW['mlog10padj_err'] = - np.log10(self.dGW['cilow_gene_padj']) \
                                     - self.dGW['mlog10padj']
        for gid in self.dGW[self.id_type].unique():
            df = self.dGW[self.dGW[self.id_type]==gid]
            gsymbol = df['hgnc_symbol'].unique()[0]
            self._barplot(df, gid, gsymbol, dom='GO')
            #Barplots separated by go domain
#             for go_dom in self.go_domains:
#                 self._barplot(df[df['go_domain']==go_dom],gid, gsymbol,dom=go_dom)
    
    
    def _get_scatter_data(self):
        """
        Data processing function to for scatter plots.
        """
        scd = dict()
        if not self.std_id
            scat_cols = [self.id_type,'hgnc_symbol','hgnc_id','con','go_con',
                         'gene_con','rel_go','frac_rel_go']
        else:
            scat_cols = [self.id_type,'hgnc_symbol','con','go_con',
                         'gene_con','rel_go','frac_rel_go']
        for gid in self.dGW[self.id_type].unique():
            df = self.dGW[ self.dGW[self.id_type] == gid ]
            gname = df['hgnc_symbol'].unique()[0]
            con = df['ncon_gene'].unique()[0]
            if np.isnan(con):#no GO annotations
                go_con = np.nan
                gene_con = np.nan
                relgo = np.nan
                frac_relgo = np.nan
            else:
                go_con = len(df)
                gene_con = con - go_con
                relgo = len(df[df['global_padj'] < alpha_FDR])
                frac_relgo = float(relgo)/go_con

            if std_id_type:
                scd[gid] = [gid, gname,con,go_con, 
                            gene_con, relgo, frac_relgo]
            else:
                hid = df['hgnc_id'].unique()[0]
                scd[gid] = [gid, gname, hid, con, go_con, 
                            gene_con, relgo, frac_relgo]
        self.scatter_data = pd.DataFrame.from_dict(scd, orient='index',
                                                   columns=scat_cols)
        self.scatter_data.to_csv(os.path.join(self.path,
                                              'genewalk_scatterplots.csv'),
                                 index=False)
        
    
    def _barplot(dplot,gene_id,gsymbol,dom):
        """
        Bar plot figure generating function.
        """
        if len(dplot)>0:
            sns.set(style="whitegrid", color_codes=True)
            f, ax = plt.subplots(figsize=(2,0.25*len(dplot))) 
            ax = sns.barplot(x='mlog10padj', y='go_name',
                             xerr=dplot['mlog10padj_err'],
                             data=dplot,color="b",
                             error_kw=dict(elinewidth=1,
                                           ecolor=[0.7,0.7,0.7],capsize=1))
            plt.axvline(x=-np.log10(self.alpha_FDR), color='r', linestyle='--')
    #         sns.set(font_scale=3)
            font_sz=12
            max_x_data = max(dplot['mlog10padj']+dplot['mlog10padj_err'])+0.3
            plt.xticks(range(round(max_x_data)+1), size=font_sz)
            plt.xlim(0, round(max_x_data,1))
            plt.xlabel('-log_{10}('+self.stat+')', size=font_sz)
            ax.yaxis.set_label_position('right')
            plt.ylabel(dom+' annotations', size=font_sz)
            plt.title(gsymbol+' '+gene_id, size=font_sz)
            filename = 'barplot_'+gsymbol+' '+gene_id+'_x_mlog10'+ \
                        self.stat+'_y_'+dom+'.'
#             plt.savefig(os.path.join(self.path,filename+'pdf'),
#                         bbox_inches="tight",transparent=True)
            plt.savefig(os.path.join(self.path,filename+'png'),
                        bbox_inches="tight",transparent=True)
        else:
            logger.warning('No results for gene id %s: \
                            could not produce barplot' % gene_id)