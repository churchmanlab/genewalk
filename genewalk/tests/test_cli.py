import os
import glob
import shutil
from genewalk.cli import run_main, default_base_folder

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_BASE_FOLDER = os.path.join(default_base_folder, '.test')


class ArgparseMock:
    def __init__(self, project, genes, id_type,
                 stage='all', base_folder=TEST_BASE_FOLDER,
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


def _place_files():
    test_resource_files = glob.glob(os.path.join(HERE, 'resources', '*'))
    test_resource_folder = \
        os.path.join(TEST_BASE_FOLDER, 'resources')
    os.mkdir(test_resource_folder)
    for test_file in test_resource_files:
        shutil.copy(test_file,
                    os.path.join(test_resource_folder,
                                 os.path.basename(test_file)))


def test_default():
    project_name = 'test1'
    gene_list = os.path.join(HERE, 'resources', 'hgnc_symbols.txt')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol')

    _place_files()

    run_main(args)

    assert os.path.exists(TEST_BASE_FOLDER)
    assert os.path.exists(os.path.join(TEST_BASE_FOLDER, project_name,
                                       'genewalk_results.csv'))
