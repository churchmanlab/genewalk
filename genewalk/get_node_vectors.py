#!/usr/bin/env python

import networkx as nx
import indra
import pickle as pkl
import copy
import time
import argparse
from genewalk.nx_mg_assembler import Nx_MG_Assembler
from genewalk.deepwalk import DeepWalk


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description='Choose a path where GeneWalk files are generated (default: ~/genewalk/ ).')
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--stmts', default='data/JQ1_HGNCidForINDRA_stmts.pkl')
    parser.add_argument('--fplx', default='data/JQ1_HGNCidForINDRA_fplx.txt')
    parser.add_argument('--path_GO', default='~/genewalk/GO/')
    parser.add_argument('--Nreps', default=10)
    args = parser.parse_args()

    # Open pickled statements
    print('loading', args.stmts)
    with open(args.path+args.stmts, 'rb') as f:
        stmts=pkl.load(f)
    
    print('assembling network')
    MG=Nx_MG_Assembler(stmts,args.path_GO)
    del(stmts)
    
    print('adding genes nodes from INDRA stmts')
    MG.MG_from_INDRA()
    MG.add_FPLXannotations(args.fplx)
    print('number of gene nodes',nx.number_of_nodes(MG.graph))
    print('adding GO nodes')
    MG.add_GOannotations()
    MG.add_GOontology()
    print('total number of nodes in network',nx.number_of_nodes(MG.graph))
    
    #pickle the network
    filename='GeneWalk_MG.pkl'
    MGA=copy.deepcopy(MG.graph)
    with open(args.path+filename, 'wb') as f:
        pkl.dump(MGA,f,protocol=pkl.HIGHEST_PROTOCOL)
    del(MG)
    
    for rep in range(1,args.Nreps+1):
        print(rep,'/',args.Nreps)
        DW=DeepWalk(graph=MGA)

        print('generate random walks')
        start = time.time()
        DW.get_walks()
        end = time.time()
        print('DW.get_walks done ',end - start)#in sec

        print('generate node vectors')
        start = time.time()
        DW.word2vec()
        end = time.time()
        print('DW.word2vec done',end - start)#in sec

        ### Pickle the node vectors (embeddings) and DW object
        nv = copy.deepcopy(DW.model.wv)
        filename='GeneWalk_DW_nv_'+str(rep)+'.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(nv,f)
        filename='GeneWalk_DW_'+str(rep)+'.pkl'
        with open(args.path+filename, 'wb') as f:
            pkl.dump(DW,f)
