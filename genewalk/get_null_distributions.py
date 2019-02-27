#!/usr/bin/env python

import networkx as nx
import pickle as pkl
import copy
import time
import random
import argparse
from genewalk.genewalk.deepwalk import DeepWalk

def Get_Rand_Graph(mg,seed):
    """argument: graph = the original nx (multi)graph. 
    -Determines degree of multigraph for each node
    -Get a random multigraph. 
    """
    d_seq=[]
    g_view=nx.nodes(mg)#this is not randomized: order is same
    for n in g_view:
        d_seq.append(mg.degree(n))
    d_seq=sorted(d_seq,reverse=True)
    R_graph=nx.configuration_model(d_seq,seed=seed)#creates random multigraph with same degree sequence
    #the node labels are numbers which gives problems in word2vec so adjust 0 to n0
    mapping=dict()
    for n in R_graph.nodes():
        mapping[n]='n'+str(n)
    R_graph=nx.relabel_nodes(R_graph,mapping, copy=False)
    return R_graph



if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description='Choose a path where GeneWalk files are generated (default: ~/genewalk/ ).')
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--mg', default='GeneWalk_MG.pkl')
    parser.add_argument('--Nreps', default=15)
    args = parser.parse_args()

    #load multigraph
    with open(args.path+args.mg, 'rb') as f:
        MGA=pkl.load(f)

    for rep in range(1,args.Nreps+1):
        #generate random graph
        random.seed(a=None)
        gseed=random.randint(0, 1e12)
        print('gseed',gseed)
        rg=Get_Rand_Graph(MGA,gseed)
        
        #initialize DW object
        DW=DeepWalk(graph=rg)
        del(rg)

        print('generate random walks', rep)
        start = time.time()
        DW.get_walks()
        end = time.time()
        print('DW.get_walks done',end - start)#in sec

        #pickling to save walks
        filename='GeneWalk_DW_rand_'+rep+'.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(DW,f)

        print('generate node vectors', rep)
        start = time.time()
        DW.word2vec()
        end = time.time()
        print('DW.word2vec ',end - start)#in sec

        ##Pickle the node vectors (embeddings) and DeepWalk object
        nv = copy.deepcopy(DW.model.wv)
        filename='GeneWalk_DW_nv_rand_'+rep+'.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(nv,f)

        filename='GeneWalk_DW_rand_'+rep+'.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(DW,f)
