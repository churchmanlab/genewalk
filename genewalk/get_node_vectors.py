"""
This module runs DeepWalk on a given network.
"""
import logging
from genewalk.deepwalk import DeepWalk

logger = logging.getLogger('genewalk.get_node_vectors')


def run_repeat(MGA):
    """Run a single repeat of DeepWalk on the given GeneWalk network.

    This function creates two pickle files, one for the node vectors,
    and one for the DeepWalk instance for reproducibility.

    Parameters
    ----------
    MGA : networkx.MultiGraph
        The GeneWalk network to run DeepWalk on.
    """
    DW = DeepWalk(graph=MGA)
    DW.get_walks()
    DW.word2vec()
    return DW
