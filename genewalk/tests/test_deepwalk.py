from networkx import MultiGraph
from genewalk.deepwalk import DeepWalk


def test_run_walks():
    mg = MultiGraph()
    mg.add_edge('a', 'b')
    mg.add_edge('b', 'c')
    dw = DeepWalk(mg, 10, 100)
    dw.get_walks(workers=1)
    # Number of neighbors for all nodes together is 4, times niter: 400
    assert len(dw.walks) == 400, len(dw.walks)
    assert dw.walks[0][0] == 'a'

    dw.get_walks(workers=2)
    # Number of neighbors for all nodes together is 4, times niter: 400
    assert len(dw.walks) == 400, len(dw.walks)
    assert dw.walks[0][0] == 'a'
