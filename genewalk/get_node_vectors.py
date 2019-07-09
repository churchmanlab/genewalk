"""
This module runs DeepWalk on a given network.
"""
import os
import copy
import time
import logging
import argparse
import functools
import pickle as pkl
import networkx as nx
from multiprocessing import Pool
from genewalk.deepwalk import DeepWalk

logger = logging.getLogger('genewalk.get_node_vectors')


def run_repeat(rep, MGA):
    """Run a single repeat of DeepWalk on the given GeneWalk network.

    This function creates two pickle files, one for the node vectors,
    and one for the DeepWalk instance for reproducibility.

    Parameters
    ----------
    rep : int
        The ID of the repeat to run.
    MGA : networkx.MultiGraph
        The GeneWalk network to run DeepWalk on.
    """
    # TODO: create file logger handle inside the project folder for each
    #  repetition separately
    logger.info('%s/%s' % (rep, args.Nreps))
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

    # Pickle the node vectors (embeddings) and DW object
    nv = copy.deepcopy(DW.model.wv)
    filename = 'GeneWalk_DW_nv_%d.pkl' % rep
    with open(os.path.join(args.path, filename), 'wb') as f:
        pkl.dump(nv, f)
    filename = 'GeneWalk_DW_%d.pkl' % rep
    with open(os.path.join(args.path, filename), 'wb') as f:
        pkl.dump(DW, f)
