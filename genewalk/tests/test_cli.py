import os
import pandas
import logging
from nose.tools import raises
from genewalk.cli import run_main
from .util import place_resource_files, TEST_RESOURCES, TEST_BASE_FOLDER

logger = logging.getLogger(__name__)


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


def test_default():
    project_name = 'test1'
    gene_list = os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol')
    place_resource_files()
    run_main(args)

    assert os.path.exists(TEST_BASE_FOLDER)
    result_csv = os.path.join(TEST_BASE_FOLDER, project_name,
                              'genewalk_results.csv')
    assert os.path.exists(result_csv)
    df = pandas.read_csv(result_csv)
    assert 'MAP2K2' in set(df['hgnc_symbol']), df['hgnc_symbol']


def test_sif():
    project_name = 'test_sif'
    gene_list = os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt')
    sif = os.path.join(TEST_RESOURCES, 'test_sif.sif')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol',
                        network_source='sif', network_file=sif)
    place_resource_files()
    run_main(args)
    assert os.path.exists(TEST_BASE_FOLDER)
    result_csv = os.path.join(TEST_BASE_FOLDER, project_name,
                              'genewalk_results.csv')
    assert os.path.exists(result_csv)
    df = pandas.read_csv(result_csv)
    assert 'MAP2K2' in set(df['hgnc_symbol']), df['hgnc_symbol']
    assert 'GO:0005515' in set(df['go_id'])
    assert 'GO:0001934' in set(df['go_id'])
    assert 'biological process' in set(df['go_domain'])


def test_edge_list():
    project_name = 'test_edge_list'
    gene_list = os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt')
    el = os.path.join(TEST_RESOURCES, 'test_edge_list.txt')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol',
                        network_source='edge_list', network_file=el)
    place_resource_files()
    run_main(args)
    assert os.path.exists(TEST_BASE_FOLDER)
    result_csv = os.path.join(TEST_BASE_FOLDER, project_name,
                              'genewalk_results.csv')
    assert os.path.exists(result_csv)
    df = pandas.read_csv(result_csv)
    assert 'MAP2K2' in set(df['hgnc_symbol']), df['hgnc_symbol']
    assert 'GO:0005515' in set(df['go_id'])
    assert 'GO:0001934' in set(df['go_id'])
    assert 'biological process' in set(df['go_domain'])


def test_sif_annot():
    project_name = 'test_sif'
    gene_list = os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt')
    sif = os.path.join(TEST_RESOURCES, 'test_sif_annot.sif')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol',
                        network_source='sif_annot', network_file=sif)
    place_resource_files()
    run_main(args)
    assert os.path.exists(TEST_BASE_FOLDER)
    result_csv = os.path.join(TEST_BASE_FOLDER, project_name,
                              'genewalk_results.csv')
    assert os.path.exists(result_csv)
    df = pandas.read_csv(result_csv)
    assert 'MAP2K2' in set(df['hgnc_symbol']), df['hgnc_symbol']
    assert 'GO:0005515' in set(df['go_id'])
    # In this case we don't have this annotation in the SIF file
    assert 'GO:0001934' not in set(df['go_id'])
    assert 'biological process' in set(df['go_domain'])


def test_sif_full():
    project_name = 'test_sif'
    gene_list = os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt')
    sif = os.path.join(TEST_RESOURCES, 'test_sif_full.sif')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol',
                        network_source='sif_full', network_file=sif)
    place_resource_files()
    run_main(args)
    assert os.path.exists(TEST_BASE_FOLDER)
    result_csv = os.path.join(TEST_BASE_FOLDER, project_name,
                              'genewalk_results.csv')
    assert os.path.exists(result_csv)
    df = pandas.read_csv(result_csv)
    assert 'MAP2K2' in set(df['hgnc_symbol']), df['hgnc_symbol']
    assert 'GO:0005515' in set(df['go_id'])
    # In this case we don't have this annotation in the SIF file
    assert 'GO:0001934' not in set(df['go_id'])
    # It's particularly important here that we have e.g., GO domains
    # since this means that node attributes are added correctly for
    # GO terms mentioned in the SIF file
    assert 'biological process' in set(df['go_domain'])


def test_custom_genes():
    project_name = 'test_custom_genes'
    gene_list = os.path.join(TEST_RESOURCES, 'custom_gene_list.txt')
    sif = os.path.join(TEST_RESOURCES, 'test_sif_custom.sif')
    args = ArgparseMock(project_name, gene_list, 'custom',
                        network_source='sif_annot', network_file=sif)
    place_resource_files()
    run_main(args)
    assert os.path.exists(TEST_BASE_FOLDER)
    result_csv = os.path.join(TEST_BASE_FOLDER, project_name,
                              'genewalk_results.csv')
    assert os.path.exists(result_csv)
    df = pandas.read_csv(result_csv)
    assert 'CUSTOM:ABC' in set(df['custom']), df['custom']
    assert 'GO:0000186' in set(df['go_id'])
    # In this case we don't have this annotation in the SIF file
    assert 'GO:0001934' not in set(df['go_id'])
    assert 'biological process' in set(df['go_domain'])


@raises(ValueError)
def test_missing_network_file():
    project_name = 'test_missing_network_file'
    gene_list = os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt')
    args = ArgparseMock(project_name, gene_list, 'hgnc_symbol',
                        network_source='sif')
    place_resource_files()
    run_main(args)
