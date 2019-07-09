"""
This module runs DeepWalk on a given network.
"""
import time
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
    # TODO: create file logger handle inside the project folder for each
    #  repetition separately
    DW = DeepWalk(graph=MGA)

    logger.info('generate random walks')
    start = time.time()
    DW.get_walks()
    end = time.time()
    logger.info('DW.get_walks done %.2f' % (end - start))  # in sec

    logger.info('generate node vectors')
    start = time.time()
    DW.word2vec()
    end = time.time()
    logger.info('DW.word2vec done %.2f' % (end - start))  # in sec
    return DW
