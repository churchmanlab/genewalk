#!/usr/bin/env python
import os
import copy
import time
import logging
import argparse
import functools
import pickle as pkl
import networkx as nx
from multiprocessing import Pool
from deepwalk import DeepWalk
from nx_mg_assembler import Nx_MG_Assembler_PC, Nx_MG_Assembler_INDRA, Nx_MG_Assembler_fromUser

logger = logging.getLogger(__name__)

def run_repeat(rep, MGA):
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



if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description='Choose a folder where files will be generated (default: ~/genewalk/ ). '
        'Define filename of list with genes of interest (default: gene_list.csv). '
        'Decide which data_source is used: Pathway Commons (default: PC), indra, '
        'or a user-provided network from file (fromUser). '
        'Set mouse_genes to True (default: False) if the gene_list contains MGI identifiers '
        'instead of human genes. '
        'A GeneWalk Network is then assembled and network representation learning performed.')
    parser.add_argument('--path', default='~/genewalk/')###############This needs work: talk with Ben what is best set up
    parser.add_argument('--genes', default='gene_list.txt')
    parser.add_argument('--data_source', default='PC')
    parser.add_argument('--mouse_genes', default=False)
    parser.add_argument('--nproc', default=1, type=int)
    parser.add_argument('--Nreps', default=10, type=int)
    args = parser.parse_args()
    
    logger.info('initializing network')
    if args.data_source == 'PC':
        MG=Nx_MG_Assembler_PC(os.path.join(args.path, args.genes),mouse_genes=args.mouse_genes)
        
        logger.info('adding gene nodes from Pathway Commons')
        MG.MG_from_PC()
        logger.info('number of PC originating nodes %d' % nx.number_of_nodes(MG.graph))
        
        logger.info('adding GO nodes')
        MG.add_GOannotations()
        MG.add_GOontology()
        
    elif args.data_source == 'indra':
        logger.info('Currently this option is not yet available for any user-provided gene list, but it will become available for public use in the future. Choose PC as data_source instead to run GeneWalk with a user-provided gene list. Now, we proceed to demonstrate the indra option by running GeneWalk on the JQ1 study described in Ietswaart et al.' )
        
        # Open pickled statements
        fstmts='INDRA_stmts.pkl'
        logger.info('loading %s' % fstmts)
        with open(os.path.join(args.path, fstmts), 'rb') as f:
            stmts = pkl.load(f)

        MG=Nx_MG_Assembler_INDRA(stmts)
        del stmts
        
        logger.info('adding nodes from INDRA stmts')
        MG.MG_from_INDRA()
        
        ffplx='INDRA_fplx.txt'
        MG.add_FPLXannotations(os.path.join(args.path, ffplx))
        
        logger.info('number of INDRA originating nodes %d' % nx.number_of_nodes(MG.graph))
        
        logger.info('adding GO nodes')
        MG.add_GOannotations()
        MG.add_GOontology()
        
    elif args.data_source == 'fromUser':
        logger.info('loading user-provided GeneWalk Network from %s' % args.genes)
        MG=Nx_MG_Assembler_fromUser(os.path.join(args.path, args.genes))
    
    else: 
        logger.info('Please specify data_source flag as PC, indra or fromUser')

    logger.info('total number of nodes in network %d' % nx.number_of_nodes(MG.graph))

    #pickle the network
    fmg='GeneWalk_MG.pkl'
    MGA=copy.deepcopy(MG.graph)
    with open(os.path.join(args.path, fmg), 'wb') as f:
        pkl.dump(MGA,f,protocol=pkl.HIGHEST_PROTOCOL)
    del MG 

    pool = Pool(args.nproc) if args.nproc > 1 else None
    if pool:
        run_repeat_wrapper = functools.partial(run_repeat, MGA=MGA)
        pool.map(run_repeat_wrapper, range(1, args.Nreps + 1))
    else:
        for rep in range(1, args.Nreps + 1):
            run_repeat(rep, MGA)
