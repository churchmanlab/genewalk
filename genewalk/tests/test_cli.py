import os
import shutil
from genewalk.cli import run_main, default_base_folder

HERE = os.path.dirname(os.path.abspath(__file__))
go_obo_path = os.path.join(default_base_folder, 'go.obo')
go_test_obo_path = os.path.join(HERE, 'go_test.obo')


class ArgparseMock:
    def __init__(self, project, genes, id_type,
                 stage='all', base_folder=default_base_folder,
                 network_source='pc', network_file=None, nproc=1,
                 nreps_graph=3, nreps_null=3, alpha_fdr=1.0,
                 dim_rep=8, save_dw=False, random_seed=None):
        self.project = project
        self.genes = genes
        self.id_type = id_type
        self.stage = stage
        self.base_folder = base_folder
        self.network_source = network_source
        self.network_file = network_file
        self.nproc = nproc
        self.nreps_graph = nreps_graph
        self.nreps_null = nreps_null
        self.alpha_fdr = alpha_fdr
        self.dim_rep = dim_rep
        self.save_dw = save_dw
        self.random_seed = random_seed


def test_default():
    project_name = 'test1'
    gene_list = os.path.join(HERE, 'hgnc_symbols.txt')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol')

    shutil.copy(go_test_obo_path, go_obo_path)

    run_main(args)

    assert os.path.exists(default_base_folder)
    assert os.path.exists(os.path.join(default_base_folder, project_name,
                                       'genewalk_results.csv'))
