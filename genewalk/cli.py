import gc
import os
import sys
import copy
import pickle
import random
import logging
import argparse
import numpy as np
from genewalk import __version__
from genewalk.nx_mg_assembler import load_network
from genewalk.gene_lists import read_gene_list
from genewalk.deepwalk import run_walks
from genewalk.null_distributions import get_rand_graph, \
    get_null_distributions
from genewalk.perform_statistics import GeneWalk
from genewalk import logger as root_logger, default_logger_format, \
    default_date_format
from genewalk.resources import ResourceManager

logger = logging.getLogger('genewalk.cli')

default_base_folder = os.path.join(os.path.expanduser('~/'), 'genewalk')


def create_project_folder(base_folder, project):
    project_folder = os.path.join(base_folder, project)
    logger.info('Creating project folder at %s' % project_folder)
    if not os.path.exists(project_folder):
        os.makedirs(project_folder)
    return project_folder


def save_pickle(obj, project_folder, prefix):
    fname = os.path.join(project_folder, '%s.pkl' % prefix)
    logger.info('Saving into %s...' % fname)
    with open(fname, 'wb') as fh:
        pickle.dump(obj, fh)


def load_pickle(project_folder, prefix):
    fname = os.path.join(project_folder, '%s.pkl' % prefix)
    logger.info('Loading %s...' % fname)
    with open(fname, 'rb') as fh:
        return pickle.load(fh)


def main():
    parser = argparse.ArgumentParser(
        description='Run GeneWalk on a list of genes provided in a text '
                    'file.')
    parser.add_argument('--version', action='version',
                        version='GeneWalk %s' % __version__,
                        help='Print the version of GeneWalk and exit.')
    parser.add_argument('--project', help='A name for the project which '
                                          'determines the folder within the '
                                          'base folder in which the '
                                          'intermediate and final results '
                                          'are written. Must contain only '
                                          'characters that are valid in '
                                          'folder names.',
                        required=True)
    parser.add_argument('--genes', help='Path to a text file with a list of '
                                        'genes of interest, for example'
                                        'differentially expressed genes. '
                                        'The type of gene identifiers used in '
                                        'the text file are provided in the '
                                        'id_type argument.',
                        required=True)
    parser.add_argument('--id_type',
                        help='The type of gene IDs provided in the text file '
                             'in the genes argument. Possible values are: '
                             'hgnc_symbol, hgnc_id, ensembl_id, and mgi_id.',
                        choices=['hgnc_symbol', 'hgnc_id',
                                 'ensembl_id', 'mgi_id', 'entrez_human',
                                 'entrez_mouse'],
                        required=True)
    parser.add_argument('--stage', default='all',
                        help='The stage of processing to run. Default: '
                             '%(default)s',
                        choices=['all', 'node_vectors', 'null_distribution',
                                 'statistics'])
    parser.add_argument('--base_folder', default=default_base_folder,
                        help='The base folder used to store GeneWalk '
                             'temporary and result files for a given project.'
                             ' Default: %(default)s')
    parser.add_argument('--network_source', default='pc',
                        help='The source of the network to be used.'
                             'Possible values are: pc, indra, edge_list, and '
                             'sif. In case of indra, edge_list, and sif, '
                             'the network_file argument must be specified.'
                             ' Default: %(default)s',
                        choices=['pc', 'indra', 'edge_list', 'sif'])
    parser.add_argument('--network_file', default=None,
                        help='If network_source is indra, this argument '
                             'points to a Python pickle file in which a list '
                             'of INDRA Statements constituting the network '
                             'is contained. In case network_source is '
                             'edge_list or sif, '
                             'the network_file argument points to a text file '
                             'representing the network.')
    parser.add_argument('--nproc', default=1, type=int,
                        help='The number of processors to use in a '
                             'multiprocessing environment. Default: '
                             '%(default)s')
    parser.add_argument('--nreps_graph', default=3, type=int,
                        help='The number of repeats to run when calculating '
                             'node vectors on the GeneWalk graph. '
                             'Default: %(default)s')
    parser.add_argument('--nreps_null', default=3, type=int,
                        help='The number of repeats to run when calculating '
                             'node vectors on the random network graphs '
                             'for constructing the null distribution. '
                             'Default: %(default)s')
    parser.add_argument('--alpha_fdr', default=1, type=float,
                        help='The false discovery rate to use when '
                             'outputting the final statistics table. '
                             'If 1 (default), all similarities are output, '
                             'otherwise only the ones whose false discovery '
                             'rate are below this parameter are included. '
                             'Default: %(default)s')
    parser.add_argument('--save_dw', default=False, type=bool,
                        help='If True, the full DeepWalk object for each '
                             'repeat is saved in the project folder. This can '
                             'be useful for debugging but the files are '
                             'typically very large. Default: %(default)s')
    parser.add_argument('--random_seed', default=None, type=int,
                        help='If provided, the random number generator is '
                             'seeded with the given value. This should only '
                             'be used if the goal is to deterministically '
                             'reproduce a prior result obtained with the same '
                             'random seed.')
    args = parser.parse_args()

    # Now we run the relevant stage of processing
    project_folder = create_project_folder(args.base_folder, args.project)

    # Add a logger specific to the project and processing stage
    log_file = os.path.join(project_folder, 'genewalk_%s.log' % args.stage)
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    project_log_handler = logging.FileHandler(log_file)
    project_log_handler.setFormatter(formatter)
    root_logger.addHandler(project_log_handler)

    if args.random_seed:
        logger.info('Running with random seed %d' % args.random_seed)
        random.seed(a=int(args.random_seed))

    # Make sure we have all the resource files
    rm = ResourceManager(base_folder=args.base_folder)
    rm.download_all()

    if args.stage in ('all', 'node_vectors'):
        genes = read_gene_list(args.genes, args.id_type, rm)
        save_pickle(genes, project_folder, 'genes')
        MG = load_network(args.network_source, args.network_file, genes,
                          resource_manager=rm)
        save_pickle(MG.graph, project_folder, 'multi_graph')
        for i in range(args.nreps_graph):
            logger.info('%s/%s' % (i + 1, args.nreps_graph))
            DW = run_walks(MG.graph, workers=args.nproc)

            # Pickle the node vectors (embeddings) and DW object
            if args.save_dw:
                save_pickle(DW, project_folder, 'deepwalk_%d' % (i + 1))
            nv = copy.deepcopy(DW.model.wv)
            save_pickle(nv, project_folder,
                        'deepwalk_node_vectors_%d' % (i + 1))

            # Delete the DeepWalk object to clear memory
            del DW, nv
            gc.collect()

    if args.stage in ('all', 'null_distribution'):
        MG = load_pickle(project_folder, 'multi_graph')
        srd = []
        for i in range(args.nreps_null):
            logger.info('%s/%s' % (i + 1, args.nreps_null))
            RG = get_rand_graph(MG)
            DW = run_walks(RG, workers=args.nproc)

            # Pickle the node vectors (embeddings) and DW object
            if args.save_dw:
                save_pickle(DW, project_folder, 'deepwalk_rand_%d' % (i + 1))
            nv = copy.deepcopy(DW.model.wv)
            save_pickle(nv, project_folder, 'deepwalk_node_vectors_rand_%d'
                                            % (i + 1))
            # Delete the DeepWalk object to clear memory
            del DW
            gc.collect()

            # Calculate the null distributions
            srd += get_null_distributions(RG, nv)
            del nv
            gc.collect()
        srd = np.asarray(sorted(srd))
        save_pickle(srd, project_folder, 'genewalk_rand_simdists')

    if args.stage in ('all', 'statistics'):
        MG = load_pickle(project_folder, 'multi_graph')
        genes = load_pickle(project_folder, 'genes')
        nvs = [load_pickle(project_folder,
                           'deepwalk_node_vectors_%d' % (i + 1))
               for i in range(args.nreps_graph)]
        null_dist = load_pickle(project_folder, 'genewalk_rand_simdists')
        GW = GeneWalk(MG, genes, nvs, null_dist)
        df = GW.generate_output(alpha_fdr=args.alpha_fdr,
                                base_id_type=args.id_type)
        fname = os.path.join(project_folder, 'genewalk_results.csv')
        logger.info('Saving final results into %s' % fname)
        df.to_csv(fname, index=False, float_format='%.3e')


if __name__ == '__main__':
    main()
