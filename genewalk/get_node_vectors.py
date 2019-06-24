#!/usr/bin/env python

import copy
import time
import argparse
import pickle as pkl
import networkx as nx
from genewalk.nx_mg_assembler import Nx_MG_Assembler_PC, Nx_MG_Assembler_INDRA, Nx_MG_Assembler_fromUser 
from genewalk.deepwalk import DeepWalk


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description='Choose a folder where files will be generated (default: ~/genewalk/ ). \
        Define filename of list with genes of interest (default: gene_list.txt). \
        Decide which data_source is used: Pathway Commons (default: PC), \
        indra, or a user-provided network from file (fromUser). \
        Set mouse_genes to True (default: False) if the gene_list contains MGI identifiers instead of human genes. \ 
        A GeneWalk Network is then assembled and network representation learning performed.')
    parser.add_argument('--path', default='~/genewalk/')###############This needs work: talk with Ben what is best set up
    parser.add_argument('--genes', default='gene_list.txt')
    parser.add_argument('--data_source', default='PC')
    parser.add_argument('--mouse_genes', default=False)
    parser.add_argument('--Nreps', default=10)
    args = parser.parse_args()
    self.path=path
    
    print('initializing network')
    if args.data_source == 'PC':
        MG=Nx_MG_Assembler_PC(args.genes)
        
        print('adding gene nodes from Pathway Commons')
        MG.MG_from_PC()
        print('number of PC originating nodes',nx.number_of_nodes(MG.graph))
        
        print('adding GO nodes')
        MG.add_GOannotations()
        MG.add_GOontology()
        
    elif args.data_source == 'indra':
        print('Currently this option is not yet available for any user-provided gene list, but it will \
                become available for public use in the future. Choose PC as data_source instead \
                to run GeneWalk with a user-provided gene list. \
                Now, we proceed to demonstrate the indra option by running GeneWalk on the JQ1 study \ 
                described in Ietswaart et al.' )
        fstmts='./genewalk/JQ1_HGNCidForINDRA_stmts.pkl'
        print('loading', fstmts)
        with open(fstmts, 'rb') as f:
            stmts=pkl.load(f)

        MG=Nx_MG_Assembler_INDRA(stmts)
        del(stmts)
        
        print('adding nodes from INDRA stmts')
        MG.MG_from_INDRA()
        
        ffplx='./genewalk/JQ1_HGNCidForINDRA_fplx.txt'
        MG.add_FPLXannotations(ffplx)
        
        print('number of INDRA originating nodes',nx.number_of_nodes(MG.graph))
        
        print('adding GO nodes')
        MG.add_GOannotations()
        MG.add_GOontology()
        
    elif args.data_source == 'fromUser':
        print('loading user-provided GeneWalk Network from ', args.genes)
        MG=Nx_MG_Assembler_fromUser(args.genes)
    
    else: 
        print('Please specify data_source flag as PC, indra or fromUser')

    
    print('total number of nodes in network',nx.number_of_nodes(MG.graph))

    #pickle the network
    filename='GeneWalk_MG.pkl'
    MGA=copy.deepcopy(MG.graph)
    with open(self.path+filename, 'wb') as f:
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
        with open(self.path+filename, 'wb') as f:
            pkl.dump(nv,f)
        filename='GeneWalk_DW_'+str(rep)+'.pkl'
        with open(self.path+filename, 'wb') as f:
            pkl.dump(DW,f)
