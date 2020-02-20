"""This module implements functions related to the construction of a null
distribution for GeneWalk networks."""
import logging
import networkx as nx

logger = logging.getLogger('genewalk.get_null_distributions')


def get_rand_graph(mg):
    """Return a random graph with the same degree distribution as the input.

    Parameters
    ----------
    mg : networkx.MultiGraph
        An input graph based on which a random graph is generated.

    Returns
    -------
    networkx.MultiGraph
        A random graph whose degree distribution matches that of the output.
    """
    # this is not randomized: order is same
    d_seq = sorted([mg.degree(n) for n in mg.nodes()], reverse=True)
    # creates random multigraph with same degree sequence
    rg = nx.configuration_model(d_seq)
    # the node labels are numbers which gives problems in word2vec
    # so adjust 0 to n0
    mapping = {n: 'n%s' % n for n in rg.nodes()}
    rg = nx.relabel_nodes(rg, mapping, copy=False)
    return rg


def get_null_distributions(rg, nv):
    """Return a distribution with similarity values between (random) node
    vectors originating from the input randomized graph.
    """
    srd = []
    # Generate null distributions from random node vectors
    for node in nx.nodes(rg):
        connections = list(rg[node])
        # only get sims for non self connections
        if node in connections:
            connections.remove(node)
        n_con_source = len(connections)
        if n_con_source > 0:
            sim_dist = list(1-nv.distances(node, other_words=connections))
            srd += sim_dist
    return srd
