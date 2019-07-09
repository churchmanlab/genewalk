"""This module constructs a set of random graphs, for a given number of
repeats, where the degree distribution matches that of an input GeneWalk
network. For each random graph, DeepWalk is run, and the null distribution
of cosine similarities is constructed from all DeepWalk outputs."""
import random
import logging
import numpy as np
import networkx as nx

logger = logging.getLogger('genewalk.get_null_distributions')


def get_rand_graph(mg):
    random.seed(a=None)
    gseed = random.randint(0, 1e12)
    logger.info('gseed: %s' % gseed)
    # this is not randomized: order is same
    d_seq = sorted([mg.degree(n) for n in mg.nodes()], reverse=True)
    # creates random multigraph with same degree sequence
    rg = nx.configuration_model(d_seq, seed=gseed)
    # the node labels are numbers which gives problems in word2vec
    # so adjust 0 to n0
    mapping = {n: 'n%s' % n for n in rg.nodes()}
    rg = nx.relabel_nodes(rg, mapping, copy=False)
    return rg


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
        if N_con_source > 0:
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


def get_srd(srs):
    # similarity (random) null distributions
    # key structure: of srd: 'd'+str(#connections of a node)
    srd = {}
    for sr in srs:
        for k, v in sr.items():
            if k in srd:
                srd[k].extend(v)
            else:
                srd[k] = v

    # To improve search speed: convert lists with similarities
    # to sorted numpy arrays.
    for skey in srd.keys():
        srd[skey] = np.asarray(sorted(srd[skey]))
    return srd

