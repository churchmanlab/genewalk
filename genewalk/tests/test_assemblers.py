from indra.statements import *
from genewalk.nx_mg_assembler import IndraNxMgAssembler


def test_indra_assembly():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    kras = Agent('KRAS', db_refs={'HGNC': '6407', 'UP': 'P01116'})
    genes = [{'HGNC': agent.db_refs['HGNC'],
              'HGNC_SYMBOL': agent.name,
              'UP': agent.db_refs['UP']} for agent in [braf, kras]]
    stmts = [Phosphorylation(kras, braf)]
    mg = IndraNxMgAssembler(genes, stmts)
    assert mg.graph
    assert mg.graph.nodes['KRAS']
    assert len(mg.graph.edges('BRAF')) > 5