import os
import pandas as pd
from genewalk.cli import create_folder, default_base_folder
from genewalk.plot import GW_Plotter


def test_generate_plots():
    project = 'plot_test'
    project_folder = create_folder(default_base_folder, project)
    dGW = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'test_results.csv'))
    figure_folder = create_folder(project_folder, 'figures')
    create_folder(figure_folder, 'barplots')
    GWp = GW_Plotter(figure_folder, dGW, 0.1)
    GWp.generate_plots()

    def file_exists(fname):
        return os.path.exists(os.path.join(figure_folder, fname))

    assert file_exists('moonlighters_x_go_con_y_frac_rel_go.png')
    assert file_exists('regulators_x_gene_con_y_frac_rel_go.png')
    assert file_exists('index.html')
    assert file_exists('barplots/barplot_BRAF_1097_x_mlog10global_padj_y_GO.png')
    assert file_exists('barplots/barplot_KRAS_6407_x_mlog10global_padj_y_GO.png')
