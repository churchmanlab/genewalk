import os
from genewalk.gene_lists import read_gene_list
from genewalk.nx_mg_assembler import UserNxMgAssembler
from genewalk.resources import ResourceManager
from .util import TEST_RESOURCES, TEST_BASE_FOLDER


rm = ResourceManager(TEST_BASE_FOLDER)
sif_genes = os.path.join(TEST_RESOURCES, 'test_sif.sif')
sif_annots = os.path.join(TEST_RESOURCES, 'test_sif_annot.sif')
sif_full = os.path.join(TEST_RESOURCES, 'test_sif_full.sif')
genes = read_gene_list(os.path.join(TEST_RESOURCES, 'hgnc_symbols.txt'),
                       id_type='hgnc_symbol', resource_manager=rm)


def test_gene_only_sif():
    mga = UserNxMgAssembler(genes, resource_manager=rm, filepath=sif_genes,
                            gwn_format='sif')
    gene_nodes = {'KRAS', 'BRAF', 'MAP2K2', 'MAPK1', 'PIK3CA', 'AKT1'}
    go_nodes = {'GO:0001934', 'GO:0005515', 'GO:0000186', 'GO:0000001',
                'GO:0032147', 'GO:0003924'}
    assert set(mga.graph.nodes()) == gene_nodes | go_nodes

    # Make sure we have GO node annotations as expected
    go_node = mga.graph.nodes['GO:0000186']
    assert go_node['GO'] == 'GO:0000186', go_node
    assert go_node['domain'] == 'biological_process', go_node
    assert go_node['name'] == 'activation of MAPKK activity', go_node

    # Make sure we have GO annotation edges
    assert ('BRAF', 'GO:0000186') in mga.graph.edges
    assert ('BRAF', 'GO:0005515') in mga.graph.edges

    # Make sure we have GO DAG edges
    assert ('GO:0000186', 'GO:0032147') in mga.graph.edges


def test_annot_sif():
    mga = UserNxMgAssembler(genes, resource_manager=rm, filepath=sif_annots,
                            gwn_format='sif_annot')
    gene_nodes = {'KRAS', 'BRAF', 'MAP2K2', 'MAPK1', 'PIK3CA', 'AKT1'}
    go_nodes = {'GO:0001934', 'GO:0005515', 'GO:0000186', 'GO:0000001',
                'GO:0032147', 'GO:0003924'}
    assert set(mga.graph.nodes()) == gene_nodes | go_nodes

    # Make sure we have GO node annotations as expected
    go_node = mga.graph.nodes['GO:0000186']
    assert go_node['GO'] == 'GO:0000186', go_node
    assert go_node['domain'] == 'biological_process', go_node
    assert go_node['name'] == 'activation of MAPKK activity', go_node

    # Make sure we have GO annotation edges
    assert ('BRAF', 'GO:0000186') in mga.graph.edges
    # This is in the test GO annotations file but not in the SIF so
    # it shouldn't be here
    assert ('BRAF', 'GO:0005515') not in mga.graph.edges

    # Make sure we have GO DAG edges
    assert ('GO:0000186', 'GO:0032147') in mga.graph.edges


def test_full_sif():
    mga = UserNxMgAssembler(genes, resource_manager=rm, filepath=sif_full,
                            gwn_format='sif_full')
    gene_nodes = {'KRAS', 'BRAF', 'MAP2K2', 'MAPK1', 'PIK3CA', 'AKT1'}
    # This time we only have the GO nodes that explicitly appear in the
    # SIF file
    go_nodes = {'GO:0005515', 'GO:0000186', 'GO:0001934'}
    assert set(mga.graph.nodes()) == gene_nodes | go_nodes

    # Make sure we have GO node annotations as expected
    go_node = mga.graph.nodes['GO:0000186']
    assert go_node['GO'] == 'GO:0000186', go_node
    assert go_node['domain'] == 'biological_process', go_node
    assert go_node['name'] == 'activation of MAPKK activity', go_node

    # Make sure we have GO annotation edges
    assert ('BRAF', 'GO:0000186') in mga.graph.edges
    # This is in the test GO annotations file but not in the SIF so
    # it shouldn't be here
    assert ('BRAF', 'GO:0005515') not in mga.graph.edges

    # Make sure we have GO DAG edges, and only ones that are in the
    # SIF, not ones in the GO ontology
    assert ('GO:0000186', 'GO:0032147') not in mga.graph.edges
    assert ('GO:0000186', 'GO:0001934') in mga.graph.edges


def test_custom_sif():
    custom_genes = read_gene_list(
        os.path.join(TEST_RESOURCES, 'custom_gene_list.txt'),
        id_type='custom', resource_manager=rm)
    sif = os.path.join(TEST_RESOURCES, 'test_sif_custom.sif')
    mga = UserNxMgAssembler(custom_genes, resource_manager=rm,
                            filepath=sif, gwn_format='sif_annot')
    # This is to make sure the network filtering works
    assert 'CUSTOM:XXXX' not in mga.graph
    assert 'CUSTOM:YYYY' not in mga.graph

    gene_nodes = {'CUSTOM:ABC', 'CUSTOM:XYZ', 'CUSTOM:ZZZ'}
    # We still have all the GO nodes that are in the test GO OBO
    # since we are in sif_annot mode here
    go_nodes = {'GO:0001934', 'GO:0005515', 'GO:0000186', 'GO:0000001',
                'GO:0032147', 'GO:0003924'}
    graph_nodes = set(mga.graph.nodes())
    assert graph_nodes == gene_nodes | go_nodes, graph_nodes

    # Make sure we have GO node annotations as expected
    go_node = mga.graph.nodes['GO:0000186']
    assert go_node['GO'] == 'GO:0000186', go_node
    assert go_node['domain'] == 'biological_process', go_node
    assert go_node['name'] == 'activation of MAPKK activity', go_node

    # Make sure we have gene-gene edges
    assert ('CUSTOM:ABC', 'CUSTOM:XYZ') in mga.graph.edges

    # Make sure we have GO annotation edges
    assert ('CUSTOM:ABC', 'GO:0000186') in mga.graph.edges

    # Make sure we have GO DAG edges
    assert ('GO:0000186', 'GO:0032147') in mga.graph.edges
