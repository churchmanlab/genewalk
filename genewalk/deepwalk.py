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


class DeepWalk(object):
    """Perform DeepWalk (node2vec), i.e., unbiased random walk over nodes
    on an undirected networkx MultiGraph.

    Parameters
    ----------
    graph : networkx.MultiGraph
        A networkx multigraph to be used as the basis for DeepWalk.
    walk_length : Optional[int]
        The length of each random walk on the graph. Default: 10
    niter : Optional[int]
        The number of iterations for each node to run (this is multiplied by
        the number of neighbors of the node when determining the overall number
        of walks to start from a given node). Default: 100

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
        """Generates walks (sentences) sampled by an (unbiased) random walk
        over the networkx MultiGraph.

        Parameters
        ----------
        workers : Optional[int]
            The number of workers to use when running random walks. If greater
            than 1, multiprocessing is used to speed up random walk generation.
            Default: 1
        """
        logger.info('Running random walks...')
        self.walks = []
        start = time.time()
        nodes = nx.nodes(self.graph)
        # In case we don't parallelize
        if workers == 1:
            for count, node in enumerate(nodes):
                walks = run_walks_for_node(node, self.graph, self.niter,
                                           self.wl)
                self.walks.extend(walks)
                if (count + 1) % 100 == 0:
                    logger.info('Walks for %d/%d nodes complete in %.2fs' %
                                (count + 1, len(nodes), time.time() - start))
        # In case we parallelize
        else:
            pool = multiprocessing.Pool(workers)
            start_nodes = get_start_nodes(self.graph, self.niter)
            chunk_size = int(float(0.9*len(start_nodes))/workers)
            walk_fun = functools.partial(run_single_walk,
                                         graph=self.graph,
                                         length=self.wl)
            for count, res in enumerate(
                    pool.imap_unordered(walk_fun, start_nodes,
                                        chunksize=chunk_size)):
                self.walks.append(res)
                if (count + 1) % chunk_size == 0:
                    logger.info('%d/%d walks complete in %.2fs' %
                                (count + 1, len(start_nodes),
                                 time.time() - start))
            logger.debug("Closing pool...")
            pool.close()
            logger.debug("Joining pool...")
            pool.join()
            logger.debug("Pool closed and joined.")
        end = time.time()
        logger.info('Running random walks done in %.2fs' % (end - start))

    def word2vec(self, sg=1, size=8, window=1, min_count=1, negative=5,
                 workers=1, sample=0):
        """Set the model based on Word2Vec
        Source: https://radimrehurek.com/gensim/models/word2vec.html

        Note that his function sets the model attribute if the DeepWalk object
        and doesn't return a value.

        Parameters
        ----------
        sg : Optional[int] {1, 0}
            Defines the training algorithm. If 1, skip-gram is employed;
            otherwise, CBOW is used. For GeneWalk this is set to 1.
        size : Optional[int]
            Dimensionality of the node vectors. Default for GeneWalk is 8.
        window : Optional[int]
            a.k.a. context size. Maximum distance between the current and
            predicted word within a sentence. For GeneWalk this is set to 1
            to assess directly connected nodes only.
        min_count : Optional[int]
            Ignores all words with total frequency lower than this. For
            GeneWalk this is set to 0.
        negative : Optional[int]
            If > 0, negative sampling will be used, the int for negative
            specifies how many "noise words" should be drawn (usually between
            5-20). If set to 0, no negative sampling is used.
            Default for GeneWalk is 5.
        workers : Optional[int]
            Use these many worker threads to train the model (=faster training
            with multicore machines).
        sample : Optional[float]
            The threshold for configuring which higher-frequency words are
            randomly downsampled, useful range is (0, 1e-5). parameter t in eq
            5 Mikolov et al. For GeneWalk this is set to 0.
        """
        logger.info('Generating node vectors...')
        start = time.time()
        self.model = Word2Vec(sentences=self.walks, sg=sg, size=size,
                              window=window, min_count=min_count,
                              negative=negative, workers=workers,
                              sample=sample)
        end = time.time()
        logger.info('Generating node vectors done in %.2fs'
                    % (end - start))


def run_single_walk(start_node, graph, length):
    """Run a single random walk on a graph from a given start node.

    Parameters
    ----------
    graph : networks.MultiGraph
        The graph on which the random walk is to be run.
    start_node : str
        The identifier of the node from which the random walk starts.
    length : int
        The length of the random walk.

    Returns
    -------
    list of str
        A path of the given length, with each element corresponding to a node
        along the path.
    """
    path = [start_node]
    for i in range(1, length):
        start_node = random.choice(list(graph[start_node]))
        path.append(start_node)
    return path


def get_start_nodes(graph, niter):
    start_nodes = []
    for node in nx.nodes(graph):
        start_nodes += [node for _ in range(niter * len(graph[node]))]
    return start_nodes


def run_walks_for_node(node, graph, niter, walk_length):
    """Run all random walks starting from a given node.

    Parameters
    ----------
    node : str
        The identifier of the node from which the walks start.
    graph : networks.MultiGraph
        The graph on which the random walks are to be run.
    niter : int
        The number of iterations to run for gene nodes.
    walk_length : int
        The length of the walk.

    Returns
    -------
    list of list of str
        A list of random walks starting from the given node.
    """
    walks = []
    for _ in range(niter * len(graph[node])):
        walk = run_single_walk(node, graph, walk_length)
        walks.append(walk)
    return walks


def run_walks(graph, **kwargs):
    """Run random walks and get node vectors on a given graph.

    Parameters
    ----------
    graph : networkx.MultiGraph
        The graph on which random walks are going to be run and node vectors
        calculated.
    **kwargs
        Key word arguments passed as the arguments of the DeepWalk constructor,
        as well as the get_walks method and the word2vec method. See the
        DeepWalk class documentation for more information on these.

    Returns
    -------
    :py:class:`genewalk.deepwalk.DeepWalk`
        A DeepWalk instance whose walks attribute contains the list
        of random walks produced on the graph.
    """
    dw_args = {'walk_length': kwargs.pop('walk_length', default_walk_length),
               'niter': kwargs.pop('niter', default_niter)}
    DW = DeepWalk(graph, **dw_args)
    DW.get_walks(kwargs.get('workers', 1))
    DW.word2vec(**kwargs)
    return DW
