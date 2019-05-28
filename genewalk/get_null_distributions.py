import copy
import time
import random
import argparse
import numpy as np
import pickle as pkl
import networkx as nx
from genewalk.deepwalk import DeepWalk


def Get_Rand_Graph(mg, seed):
    """argument: graph = the original nx (multi)graph. 
    -Determines degree of multigraph for each node
    -Get a random multigraph. 
    """
    d_seq = []
    # this is not randomized: order is same
    g_view = nx.nodes(mg)
    for n in g_view:
        d_seq.append(mg.degree(n))
    d_seq = sorted(d_seq,reverse=True)
    # creates random multigraph with same degree sequence
    R_graph = nx.configuration_model(d_seq, seed=seed)
    # the node labels are numbers which gives problems in word2vec
    # so adjust 0 to n0
    mapping = dict()
    for n in R_graph.nodes():
        mapping[n] = 'n' + str(n)
    R_graph = nx.relabel_nodes(R_graph, mapping, copy=False)
    return R_graph


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description=('Choose a path where GeneWalk files are '
                     'generated (default: ~/genewalk/ ).'))
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--Nreps', default=15)
    args = parser.parse_args()

    # similarity (random) null distributions
    # key structure: of srd: 'd'+str(#connections of a node)
    srd = dict()

    # load multigraph
    with open(args.path+'GeneWalk_MG.pkl', 'rb') as f:
        MGA = pkl.load(f)

    for rep in range(1, args.Nreps+1):
        print(rep,'/', args.Nreps)
        # generate random graph
        random.seed(a=None)
        gseed = random.randint(0, 1e12)
        print('gseed', gseed)
        rg = Get_Rand_Graph(MGA, gseed)
        
        DW = DeepWalk(graph=rg)
        print('generate random walks', rep)
        start = time.time()
        DW.get_walks()
        end = time.time()
        print('DW.get_walks done', end - start)  # in sec
        # pickling to save walks
        filename='GeneWalk_DW_rand_' + rep + '.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(DW, f)

        print('generate node vectors', rep)
        start = time.time()
        DW.word2vec()
        end = time.time()
        print('DW.word2vec ', end - start)  # in sec
        # Pickle the node vectors (embeddings) and DeepWalk object
        nv = copy.deepcopy(DW.model.wv)
        filename='GeneWalk_DW_nv_rand_' + rep + '.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(nv,f)
        filename='GeneWalk_DW_rand_' + rep + '.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(DW,f)
        del DW
    
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

    # To improve search speed: convert lists with similarities
    # to sorted numpy arrays.
    for skey in srd.keys():
        srd[skey] = np.asarray(sorted(srd[skey]))
            
    # Pickle srd for use in genewalk.py
    filename = 'GeneWalk_DW_rand_simdists.pkl'
    print(filename)     
    with open(args.path + filename, 'wb') as f:
        pkl.dump(srd, f)
