"""
This module implements a wrapper around DeepWalk as a class. The class
contains a graph used as a basis for running the Deepwalk algorithm. It also
implements a method to run a given number of walks and save the walks as an
attribute of the instance.
"""
import time
import random
import logging
import functools
import networkx as nx
import multiprocessing
from gensim.models import Word2Vec


logger = logging.getLogger('genewalk.deepwalk')

default_walk_length = 10
default_niter = 100


def run_walk(graph, **kwargs):
    dw_args = {'walk_length': kwargs.pop('walk_length', default_walk_length),
               'niter': kwargs.pop('niter', default_niter)}
    DW = DeepWalk(graph, **dw_args)
    DW.get_walks()
    DW.word2vec(**kwargs)
    return DW


class DeepWalk(object):
    """Perform DeepWalk (node2vec), i.e., unbiased random walk over nodes
    on an undirected networkx MultiGraph.

    Parameters
    ----------
    graph : networkx.MultiGraph
        A networkx multigraph to be used as the basis for DeepWalk.
    walk_length : Optional[int]
        Default: 10
    niter : Optional[int]
        Default: 100

    Attributes
    ----------
    walks : list
        A list of walks.
    """
    def __init__(self, graph, walk_length=default_walk_length,
                 niter=default_niter):
        self.graph = graph
        self.walks = []
        self.wl = walk_length
        self.niter = niter
        self.model = None

    def get_walks(self, workers=1):
        """Generate collection of graph walks: one for each node
        (= starting point) sampled by an (unbiased) random walk over the
        networkx MultiGraph.
        """
        logger.info('Running random walks...')
        self.walks = []
        start = time.time()
        if workers == 1:
            for node in nx.nodes(self.graph):
                walks = run_walks_for_node(node, self.graph, self.niter,
                                           self.wl)
                self.walks.extend(walks)
        else:
            pool = multiprocessing.Pool(workers)
            walk_fun = functools.partial(run_walks_for_node,
                                         graph=self.graph,
                                         niter=self.niter,
                                         walk_length=self.wl)
            res = pool.map(walk_fun, nx.nodes(self.graph))
        end = time.time()
        logger.info('Running random walks done in %.2fs' % (end - start))

    # TODO: set worker size depending on the number of processors,
    #  do a benchmark on a large machine to see how much workers defined
    #  here speed up the process. If it scales well, we can do
    #  parallelization at this level only.
    def word2vec(self, sg=1, size=8, window=1, min_count=1, negative=5,
                 workers=4, sample=0):
        """Set the model based on Word2Vec
        Source: https://radimrehurek.com/gensim/models/word2vec.html

        Parameters
        ----------
        sentences : iterable of iterables
            The sentences iterable can be simply a list of lists of tokens,
            but for larger corpora, consider an iterable that streams the
            sentences directly from disk/network.
        sg : int {1, 0}
            Defines the training algorithm. If 1, skip-gram is employed;
            otherwise, CBOW is used. For GeneWalk this is set to 1.
        size : int
            Dimensionality of the node vectors. Default for GeneWalk is 8.
        window : int
            a.k.a. context size. Maximum distance between the current and
            predicted word within a sentence. For GeneWalk this is set to 1
            to assess directly connected nodes only.
        min_count : int
            Ignores all words with total frequency lower than this. For
            GeneWalk this is set to 0.
        negative : int
            If > 0, negative sampling will be used, the int for negative
            specifies how many "noise words" should be drawn (usually between
            5-20). If set to 0, no negative sampling is used.
            Default for GeneWalk is 5.
        workers : int
            Use these many worker threads to train the model (=faster training
            with multicore machines).
        sample : float
            The threshold for configuring which higher-frequency words are
            randomly downsampled, useful range is (0, 1e-5). parameter t in eq
            5 Mikolov et al. For GeneWalk this is set to 0.
        """
        logger.info('generate node vectors')
        start = time.time()
        self.model = Word2Vec(sentences=self.walks, sg=sg, size=size,
                              window=window, min_count=min_count,
                              negative=negative, workers=workers,
                              sample=sample)
        end = time.time()
        logger.info('DW.word2vec done %.2f' % (end - start))  # in sec


def run_single_walk(graph, start_node, length):
    """Generates walks (sentences) sampled by an (unbiased) random walk
    over the networkx MultiGraph: node and edge names for the sentences.

    Parameters
    ----------
    idx : int
        index of walk in self.walks that will form corpus for word2vec
    u : str
        starting node
    """
    path = [start_node]
    for i in range(1, length):
        start_node = random.choice(list(graph[start_node]))
        path.append(start_node)
    return path


def run_walks_for_node(node, graph, niter, walk_length):
    walks = []
    for _ in range(niter * len(graph[node])):
        walk = run_single_walk(graph, node, walk_length)
        walks.append(walk)
    return walks
