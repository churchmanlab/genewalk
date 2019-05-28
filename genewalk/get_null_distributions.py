import os
import copy
import time
import random
import logging
import argparse
import numpy as np
import pickle as pkl
import networkx as nx
from multiprocessing import Pool
from genewalk.deepwalk import DeepWalk


logger = logging.getLogger(__name__)


def get_rand_graph(mg, seed):
    """argument: graph = the original nx (multi)graph. 
    -Determines degree of multigraph for each node
    -Get a random multigraph. 
    """
    # this is not randomized: order is same
    d_seq = sorted([mg.degree(n) for n in mg.nodes()], reverse=True)
    # creates random multigraph with same degree sequence
    R_graph = nx.configuration_model(d_seq, seed=seed)
    # the node labels are numbers which gives problems in word2vec
    # so adjust 0 to n0
    mapping = {n: 'n%s' % n for n in R_graph.nodes()}
    R_graph = nx.relabel_nodes(R_graph, mapping, copy=False)
    return R_graph


def run_walk(rg, fname):
    DW = DeepWalk(graph=rg)
    logger.info('generate random walks')
    start = time.time()
    DW.get_walks()
    end = time.time()
    logger.info('DW.get_walks done: %s' % (end - start))  # in sec
    # pickling to save walks
    with open(os.path.join(args.path, fname), 'wb') as f:
        pkl.dump(DW, f)
    return DW


def get_word2vec(DW, w2v_file, walks_file):
    logger.info('generate node vectors')
    start = time.time()
    DW.word2vec()
    end = time.time()
    logger.info('DW.word2vec: %s' % (end - start))  # in sec
    # Pickle the node vectors (embeddings) and DeepWalk object
    nv = copy.deepcopy(DW.model.wv)
    with open(os.path.join(args.path, w2v_file), 'wb') as f:
        pkl.dump(nv, f)
    with open(os.path.join(args.path, walks_file), 'wb') as f:
        pkl.dump(DW, f)
    return DW, nv


def get_null_distributions(rg, nv):
    srd = {}
    # Generate null distributions from random node vectors
    g_view = nx.nodes(rg)
    for node in g_view:
        connections = list(rg[node])
        # only get sims for non self connections
        if node in connections:
            connections.remove(node)
        N_con_source = len(connections)
        sim_dist = list(1-nv.distances(node, other_words=connections))
        for inb in range(len(connections)):
            nb = connections[inb]
            N_con_neighbor = len(rg[nb])
            # make log2 scale categories of connectivity
            key = 'd' + str(np.floor(np.log2(min(N_con_source,
                                                 N_con_neighbor))))
            if key in srd:
                srd[key].append(sim_dist[inb])
            else:
                srd[key] = [sim_dist[inb]]
    return srd


def run_repeat(rep):
    logger.info('%s/%s' % (rep, args.Nreps))
    w2v_file = 'GeneWalk_DW_nv_rand_' + rep + '.pkl'
    walks_file = 'GeneWalk_DW_rand_' + rep + '.pkl'
    # generate random graph
    random.seed(a=None)
    gseed = random.randint(0, 1e12)
    logger.info('gseed: %s' % gseed)
    rg = get_rand_graph(MGA, gseed)
    DW = run_walk(rg, walks_file)
    DW, nv = get_word2vec(DW, w2v_file, walks_file)
    sr = get_null_distributions(rg, nv)
    return sr


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description=('Choose a path where GeneWalk files are '
                     'generated (default: ~/genewalk/ ).'))
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--Nreps', default=15)
    parser.add_argument('--nproc', default=1)
    args = parser.parse_args()

    # load multigraph
    with open(os.path.join(os.path.expanduser(args.path),
                           'GeneWalk_MG_QKI.pkl'), 'rb') as f:
        MGA = pkl.load(f)

    pool = Pool(args.nproc) if args.nproc > 1 else None

    if pool:
        srs = pool.map(run_repeat, range(1, args.Nreps + 1))
    else:
        srs = [run_repeat(rep) for rep in range(1, args.Nreps + 1)]
    # similarity (random) null distributions
    # key structure: of srd: 'd'+str(#connections of a node)
    srd = {}
    for sr in srs:
        for k, v in sr.items():
            srd[k] = v

    # To improve search speed: convert lists with similarities
    # to sorted numpy arrays.
    for skey in srd.keys():
        srd[skey] = np.asarray(sorted(srd[skey]))

    # Pickle srd for use in corecore.py
    filename = 'GeneWalk_DW_rand_simdists.pkl'
    logger.info(filename)
    with open(os.path.join(os.path.expanduser(args.path),
                           filename, 'wb')) as f:
        pkl.dump(srd, f)
