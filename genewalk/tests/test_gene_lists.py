import os
from nose.tools import raises
from genewalk.gene_lists import *
from genewalk.cli import default_base_folder
from genewalk.resources import ResourceManager
from .util import TEST_RESOURCES

rm = ResourceManager()
gm = GeneMapper(rm)


def test_map_lists():
    refs = map_hgnc_symbols(['BRAF', 'KRAS'], gm)
    assert refs[0]['HGNC'] == '1097', refs
    assert refs[0]['UP'] == 'P15056', refs
    assert refs[0]['HGNC_SYMBOL'] == 'BRAF', refs
    assert refs[1]['HGNC'] == '6407', refs
    assert refs[1]['UP'] == 'P01116', refs
    assert refs[1]['HGNC_SYMBOL'] == 'KRAS', refs

    refs = map_hgnc_ids(['1097', '6407'], gm)
    assert refs[0]['HGNC'] == '1097', refs
    assert refs[0]['UP'] == 'P15056', refs
    assert refs[0]['HGNC_SYMBOL'] == 'BRAF', refs
    assert refs[1]['HGNC'] == '6407', refs
    assert refs[1]['UP'] == 'P01116', refs
    assert refs[1]['HGNC_SYMBOL'] == 'KRAS', refs

    refs = map_mgi_ids(['MGI:892970'], gm)
    assert refs[0]['HGNC'] == '6817', refs
    assert refs[0]['HGNC_SYMBOL'] == 'MAL', refs
    assert refs[0]['UP'] == 'P21145', refs

    refs = map_ensembl_ids(['ENSG00000157764'], gm)
    assert refs[0]['HGNC'] == '1097', refs
    assert refs[0]['UP'] == 'P15056', refs
    assert refs[0]['HGNC_SYMBOL'] == 'BRAF', refs

    refs = map_ensembl_ids(['ENSG00000157764.9'], gm)
    assert refs[0]['HGNC'] == '1097', refs
    assert refs[0]['UP'] == 'P15056', refs
    assert refs[0]['HGNC_SYMBOL'] == 'BRAF', refs


def test_read_gene_list():
    with open('test_gene_list.txt', 'w') as fh:
        fh.write('HGNC:1097')
    refs = read_gene_list('test_gene_list.txt', 'hgnc_id', rm)
    assert len(refs) == 1


@raises(ValueError)
def test_read_gene_list_bad():
    with open('test_gene_list.txt', 'w') as fh:
        fh.write('HGNC:1097')
    refs = read_gene_list('test_gene_list.txt', 'ensembl_id', rm)
    assert len(refs) == 1


def test_read_gene_list_entrez_human():
    with open('test_gene_list_entrez_human.txt', 'w') as fh:
        fh.write('2597')
    refs = read_gene_list('test_gene_list_entrez_human.txt',
                          'entrez_human', rm)
    assert refs[0]['HGNC_SYMBOL'] == 'GAPDH'
    assert len(refs) == 1


def test_read_gene_list_entrez_mouse():
    rm = ResourceManager(base_folder=default_base_folder)
    with open('test_gene_list_entrez_mouse.txt', 'w') as fh:
        fh.write('14433')
    refs = read_gene_list('test_gene_list_entrez_mouse.txt',
                          'entrez_mouse', rm)
    assert len(refs) == 1
    assert refs[0]['MGI'] == '95640'
    assert refs[0]['HGNC_SYMBOL'] == 'GAPDH'


def test_read_gene_list_rgd():
    rm = ResourceManager(base_folder=default_base_folder)
    with open('test_gene_list_rgd.txt', 'w') as fh:
        fh.write('2561\n')
        fh.write('RGD:69323')
    refs = read_gene_list('test_gene_list_rgd.txt',
                          'rgd_id', rm)
    assert len(refs) == 2, refs
    assert refs[0]['RGD'] == '2561'
    assert refs[0]['HGNC_SYMBOL'] == 'ERBB2'
    assert refs[1]['RGD'] == '69323'
    assert refs[1]['HGNC_SYMBOL'] == 'ERBB3'


def test_read_custom_list():
    rm = ResourceManager(base_folder=default_base_folder)
    gene_list_file = os.path.join(TEST_RESOURCES, 'custom_gene_list.txt')
    refs = read_gene_list(gene_list_file, 'custom', rm)
    assert len(refs) == 3, refs
    assert refs[0] == {'ID': 'CUSTOM:ABC'}, refs
