import os.path
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import seaborn as sns
import plotly.express as px

logger = logging.getLogger('genewalk.plot')


class GW_Plotter(object):
    """
    genewalk plotter object to visualize genewalk output data.

    Parameters
    ----------
    path : str
        Path to output figures (Default: project_folder/figures)
    dgw : pandas.dataFrame
        GeneWalk results data frame.
    alpha_fdr : float
        FDR significance level for global analyses and barplots.
        Default: 0.1 and must be < 1 otherwise set to default.
        Lower this value to increase the stringency of the
        regulator selection procedure.
    """
    def __init__(self, path, dgw, alpha_fdr):
        self.path = path
        if alpha_fdr == 1:  # default input
            self.alpha_fdr = 0.1
        else:
            self.alpha_fdr = alpha_fdr
        self.scatter_data = pd.DataFrame()
        self.stat = 'global_padj'  # for barplots
        self.ci_stat = 'cilow_'+self.stat
        self.dGW = dgw
        self.go_domains = set(self.dGW['go_domain'])
        if self.dGW.columns[0] in {'mgi_id', 'rgd_id', 'ensembl_id',
                                   'entrez_human', 'entrez_mouse',
                                   'custom'}:
            self.std_id = False
            self.id_type = self.dGW.columns[0]
        else:
            self.std_id = True
            self.id_type = self.dGW.columns[1]

        self.name_namespace = 'custom' if self.id_type == 'custom' \
            else 'hgnc_symbol'
        self.id_namespace = 'custom' if self.id_type == 'custom' \
            else 'hgnc_id'

    def generate_plots(self):
        """Wrapper that calls scatter and bar plot generating functions."""
        reg_html = self.scatterplot_regulators()
        moonlight_html = self.scatterplot_moonlighters()
        self.barplot_goanno()
        self.make_html([reg_html, moonlight_html])

    def scatterplot_regulators(self):
        """Scatter plot with fraction of (globally) relevant GO annotations
        as a function of gene connectivity (to other genes) for all input
        genes.

        Genes with symbols listed are regulator genes.
        See genewalk_scatterplots.csv for full data and
        genewalk_regulators.csv for regulator genes of interest.
        Visualization thresholds:
        T_frac: minimal fraction of relevant GO annotations, set to 0.5
        T_gcon: minimal gene connectivity (to other
        genes), set to 75th quartile of distribution
        """
        if self.scatter_data.empty:
            self._get_scatter_data()
        plot_title = 'Regulator genes'
        xvar = 'gene_con'
        yvar = 'frac_rel_go'
        xlab = 'Connections with other genes (per gene)'
        ylab = 'Fraction of relevant GO terms (per gene)'
        xmin = 0.45
        xmax = max(self.scatter_data[xvar])*1.2
        T_gcon = self.scatter_data[xvar].quantile(q=0.75)
        T_frac = 0.5

        # seaborn static plot
        sns.set(style="whitegrid")
        fig, ax = plt.subplots(figsize=(12, 12))  # inches
        g = sns.scatterplot(x=xvar, y=yvar, hue=yvar,
                            linewidth=0, alpha=0.5,
                            sizes=(40, 400),
                            data=self.scatter_data,
                            ax=ax, legend=False)
        plt.axvline(x=T_gcon, color='grey', linestyle='--')
        plt.axhline(y=T_frac, color='grey', linestyle='--')
        font_sz = 16
        plt.xlabel(xlab, size=font_sz)
        plt.ylabel(ylab, size=font_sz)
        plt.xlim([xmin, xmax])
        plt.xticks(size=font_sz)
        plt.yticks(size=font_sz)

        regulators = []
        dreg = self.scatter_data[self.scatter_data[xvar] >= T_gcon]
        dreg = dreg[dreg[yvar] >= T_frac]
        for r in dreg.index:
            gname = dreg[self.name_namespace][r]
            regulators.append(gname)
            x_txt = dreg[xvar][r]
            y_txt = dreg[yvar][r]+np.random.normal(0, 0.002)
            g.text(x_txt, y_txt, gname, size=6, horizontalalignment='center',
                   color='black', weight='light', fontstyle='italic')
        g.set(xscale="log")
        plt.title(plot_title, size=font_sz)
        filename = 'regulators_x_' + xvar + '_y_' + yvar
        plt.savefig(os.path.join(self.path, filename + '.pdf'),
                    bbox_inches="tight", transparent=True)
        plt.savefig(os.path.join(self.path, filename + '.png'),
                    bbox_inches="tight", transparent=True)

        # plotly interactive plot
        fig = px.scatter(self.scatter_data[~self.scatter_data[yvar].isna()],
                         x=xvar, y=yvar,
                         color=yvar, size='rel_go',
                         hover_name=self.name_namespace,
                         hover_data=[self.name_namespace, self.id_type],
                         title=plot_title, labels={xvar: xlab, yvar: ylab},
                         log_x=True, range_x=[xmin, xmax])
        fig.add_shape(type='rect', x0=T_gcon, y0=T_frac, x1=xmax, y1=1,
                      fillcolor="LightSkyBlue", opacity=0.2,
                      layer="below", line_width=0)
        fig.add_shape(type='line', x0=xmin, y0=T_frac, x1=xmax, y1=T_frac,
                      line=dict(color='grey', dash='dash'))
        fig.add_shape(type='line', x0=T_gcon, y0=0, x1=T_gcon, y1=1,
                      line=dict(color='grey', dash='dash'))
        plotly_html = fig.to_html(full_html=False)
        fig.write_html(os.path.join(self.path, filename + '.html'))
        logger.info('%s plotted in %s...' % (plot_title, filename))

        df = pd.DataFrame(sorted(regulators), columns=['gw_regulators'])
        filename = 'genewalk_regulators.csv'
        df.to_csv(os.path.join(self.path, filename), index=False)
        logger.info('%s listed in %s...' % (plot_title, filename))
        return plotly_html

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
        if self.scatter_data.empty:
            self._get_scatter_data()
        plot_title = 'Moonlighting genes'
        xvar = 'go_con'
        yvar = 'frac_rel_go'
        xlab = 'Number of GO annotations (per gene)'
        ylab = 'Fraction of relevant GO terms (per gene)'
        xmin = 0.45
        xmax = max(self.scatter_data[xvar])*1.2
        T_gocon = max(30, self.scatter_data[xvar].quantile(q=0.75))
        T_frac = 0.5

        #seaborn static plot
        sns.set(style="whitegrid")
        fig, ax = plt.subplots(figsize=(12, 12))  # inches
        g = sns.scatterplot(x=xvar, y=yvar, hue=yvar,
                            linewidth=0, alpha=0.5,
                            sizes=(40, 400),
                            data=self.scatter_data,
                            ax=ax, legend=False)
        plt.axvline(x=T_gocon, color='grey', linestyle='--')
        plt.axhline(y=T_frac, color='grey', linestyle='--')
        font_sz=16
        plt.xlabel(xlab, size=font_sz)
        plt.ylabel(ylab, size=font_sz)
        plt.xlim([xmin, xmax])
        plt.xticks(size=font_sz)
        plt.yticks(size=font_sz)

        moonlighters = []
        dmoon = self.scatter_data[self.scatter_data[xvar] >= T_gocon]
        dmoon = dmoon[(dmoon[yvar] < T_frac) & (dmoon[yvar] > 0)]
        for m in dmoon.index:
            gname = dmoon[self.name_namespace][m]
            moonlighters.append(gname)
            x_txt = dmoon[xvar][m]
            y_txt = dmoon[yvar][m]
            g.text(x_txt, y_txt, gname, size=6, horizontalalignment='center',
                   color='black', weight='light',
                   fontstyle='italic')
        g.set(xscale="log")
        plt.title(plot_title, size=font_sz)
        filename = 'moonlighters_x_' + xvar + '_y_' + yvar
        plt.savefig(os.path.join(self.path, filename + '.pdf'),
                    bbox_inches="tight", transparent=True)
        plt.savefig(os.path.join(self.path, filename + '.png'),
                    bbox_inches="tight", transparent=True)

        ### plotly interactive plot
        fig = px.scatter(self.scatter_data[~self.scatter_data[yvar].isna()],
                         x=xvar, y=yvar,
                         color=yvar, size='gene_con',
                         hover_name=self.name_namespace,
                         hover_data=[self.name_namespace, self.id_type],
                         title=plot_title, labels={xvar: xlab, yvar: ylab},
                         log_x=True, range_x=[xmin, xmax])
        fig.add_shape(type='rect', x0=T_gocon, y0=0, x1=xmax, y1=T_frac,
                      fillcolor="LightSkyBlue", opacity=0.2,
                      layer="below", line_width=0)
        fig.add_shape(type='line', x0=xmin, y0=T_frac, x1=xmax, y1=T_frac,
                      line=dict(color='grey', dash='dash'))
        fig.add_shape(type='line', x0=T_gocon, y0=0, x1=T_gocon, y1=1,
                      line=dict(color='grey', dash='dash'))
        fig.write_html(os.path.join(self.path, filename + '.html'))
        plotly_html = fig.to_html(full_html=False)
        logger.info('%s plotted in %s...' % (plot_title, filename))

        df = pd.DataFrame(sorted(moonlighters), columns=['gw_moonlighter'])
        filename = 'genewalk_moonlighters.csv'
        df.to_csv(os.path.join(self.path, filename), index=False)
        logger.info('%s listed in %s...' % (plot_title, filename))
        return plotly_html

    def barplot_goanno(self):
        """Visualize statistical significances of all GO annotations in a barplot
        for a given gene of interest.
        """
        self.dGW['mlog10padj'] = -np.log10(self.dGW[self.stat])
        self.dGW['mlog10padj_err'] = - np.log10(self.dGW[self.ci_stat]) \
                                     - self.dGW['mlog10padj']
        for gid in self.dGW[self.id_type].unique():
            df = self.dGW[self.dGW[self.id_type] == gid]
            gsymbol = self.scatter_data[self.name_namespace][gid]
            self._barplot(df, gid, gsymbol, dom='GO')
            #Barplots separated by go domain
            #for go_dom in self.go_domains:
            #    self._barplot(df[df['go_domain']==go_dom], gid, gsymbol,dom=go_dom)

    def _get_scatter_data(self):
        """
        Data processing function for scatter plots.
        """
        scd = dict()
        if not self.std_id:
            if self.id_type == 'custom':
                scat_cols = [self.id_type]
            else:
                scat_cols = [self.id_type, self.name_namespace,
                             self.id_namespace]
            scat_cols += ['con', 'go_con', 'gene_con', 'rel_go',
                          'frac_rel_go']
        else:
            scat_cols = [self.id_type, self.name_namespace, 'con', 'go_con',
                         'gene_con', 'rel_go', 'frac_rel_go']
        for gid in self.dGW[self.id_type].unique():
            df = self.dGW[ self.dGW[self.id_type] == gid ]
            gname = df[self.name_namespace].unique()[0]
            con = df['ncon_gene'].unique()[0]
            if pd.isna(df['go_id'].unique()[0]):  # no GO annotations
                gocon = np.nan
                genecon = con
                relgo = 0
                fracrelgo = 0
            else:
                gocon = len(df)
                genecon = con - gocon
                relgo = min(len(df[df['global_padj'] < self.alpha_fdr]),
                            len(df[df['gene_padj'] < self.alpha_fdr]))
                fracrelgo = round(float(relgo)/gocon, 3)

            if self.std_id:
                scd[gid] = [gid, gname, con, gocon,
                            genecon, relgo, fracrelgo]
            else:
                hid = str(df[self.id_namespace].unique()[0])
                if self.id_type == 'custom':
                    scd[gid] = [gid]
                else:
                    scd[gid] = [gid, gname, hid]
                scd[gid] += [con, gocon, genecon, relgo, fracrelgo]
        self.scatter_data = pd.DataFrame.from_dict(scd, orient='index',
                                                   columns=scat_cols)
        self.scatter_data.sort_values(by=['gene_con', 'frac_rel_go', 'rel_go'],
                                      ascending=[False, False, False],
                                      inplace=True)
        filename = 'genewalk_scatterplots.csv'
        self.scatter_data.to_csv(os.path.join(self.path, filename),
                                 index=False)
        logger.info('Scatter plot data output to %s...' % filename)
        for c in ['go_con', 'gene_con']:  # for log scale plotting: 0 -> 0.5
            self.scatter_data[c].where(self.scatter_data[c] > 0, 0.5,
                                       inplace=True)

    def _barplot(self, dplot, gene_id, gsymbol, dom):
        """
        Bar plot figure generating function.
        """
        gocon = self.scatter_data['go_con'][gene_id]
        gene_id = str(gene_id)
        if gocon >= 1:
            dplot = dplot.sort_values(by=['mlog10padj', 'sim',
                                          'global_padj', 'gene_padj'],
                                      ascending=[False, False,
                                                 True, True])
            sns.set(style="whitegrid", color_codes=True)
            f, ax = plt.subplots(figsize=(2, 0.25*len(dplot)))
            ax = sns.barplot(x='mlog10padj', y='go_name',
                             xerr=dplot['mlog10padj_err'],
                             data=dplot, color="b",
                             error_kw=dict(elinewidth=1,
                                           ecolor=[0.7, 0.7, 0.7], capsize=1))
            plt.axvline(x=-np.log10(self.alpha_fdr), color='r', linestyle='--')
            font_sz = 12
            max_x_data = max(dplot['mlog10padj']+dplot['mlog10padj_err'])+0.3
            plt.xticks(range(round(max_x_data)+1), size=font_sz)
            plt.xlim(0, max(round(max_x_data, 1), 1))
            plt.xlabel('-log10('+self.stat+')', size=font_sz)
            ax.yaxis.set_label_position('right')
            plt.ylabel(dom+' annotations', size=0.8*font_sz,
                       rotation=270, labelpad=15)
            plt.title(gsymbol, size=font_sz)
            filename = 'barplot_' + gsymbol + '_' + gene_id + '_x_mlog10' + \
                self.stat + '_y_' + dom.replace(' ', '') + '.'
            #plt.savefig(os.path.join(self.path, 'barplots', filename+'pdf'),
            #            bbox_inches="tight", transparent=True)
            #pdf file size: >300kb, so memorywise not advisable to make
            #pdf barplots in case of long input gene lists.
            #Save as png (10kb) instead:
            plt.savefig(os.path.join(self.path, 'barplots', filename+'png'),
                        bbox_inches="tight", transparent=True)
            logger.info('Barplot generated for %s in %s...' %
                        (gene_id, filename))
        else:
            logger.warning('No results for gene id %s: \
                            could not produce barplot' % gene_id)

    def make_html(self, global_htmls):
        """
        Generates index.html, the html file that shows all the visualizations.
        """
        template_path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), 'results_template.html')
        with open(template_path, 'r') as fh:
            template = fh.read()

        template = template.replace('{{ GLOBAL_RESULTS }}',
                                    '\n'.join(global_htmls))

        genes = iter(self.dGW[self.id_type].unique())
        gene_results_html = ""
        while True:
            gene_results_html += '<div class="row">'
            try:
                for _ in range(4):
                    gid = next(genes)
                    gsymbol = self.scatter_data[self.name_namespace][gid]
                    img_path = os.path.join('barplots',
                               ('barplot_%s_%s_x_mlog10%s_y_GO.png' %
                                (gsymbol, gid, self.stat)))
                    # Skip any genes for which results were not generated
                    if not os.path.exists(os.path.join(self.path, img_path)):
                        continue
                    # We only make proper links for non-custom gene symbols
                    if self.id_type != 'custom':
                        identifiers_href = 'https://identifiers.org/hgnc.symbol:%s' % gsymbol
                    else:
                        identifiers_href = '#'
                    gene_results_html += """
                    <div class="col-md-3">
                        <div class="thumbnail">
                            <a href="{img_path}">
                                <img src='{img_path}' style="width:100%">
                            </a>
                            <div class="caption">
                                <a href="{identifiers_href}">{symbol} ({gid})</a>
                            </div>
                        </div>
                    </div>
                    """.format(symbol=gsymbol, gid=gid, img_path=img_path,
                               identifiers_href=identifiers_href)
            except StopIteration:
                break
            finally:
                gene_results_html += '</div>'

        template = template.replace('{{ GENE_RESULTS }}', gene_results_html)
        output_html = os.path.join(self.path, 'index.html')
        with open(output_html, 'w') as fh:
            fh.write(template)
        logger.info('index.html file generated with interactive visualizations ' 
                    'of all GeneWalk results...')
