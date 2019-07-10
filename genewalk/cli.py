import os
import copy
import pickle
import logging
import argparse
from genewalk.nx_mg_assembler import load_network
from genewalk.gene_lists import read_gene_list
from genewalk.deepwalk import run_walk
from genewalk.null_distributions import get_rand_graph, \
    get_null_distributions, get_srd
from genewalk.perform_statistics import GeneWalk


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run GeneWalk on a list of genes provided in a text file.')
    parser.add_argument('--project', help='A name for the project which '
                                          'determines the folder within the '
                                          'base folder in which the '
                                          'intermediate and final results '
                                          'are written. Must contain only '
                                          'characters that are valid in '
                                          'folder names.',
                        required=True)
    parser.add_argument('--genes', help='Path to a text file with a list of '
                                        'differentially expressed genes. The'
                                        'type of gene identifiers used in '
                                        'the text file are provided in the '
                                        'id_type argument.',
                        required=True)
    parser.add_argument('--id_type',
                        help='The type of gene IDs provided in the text file '
                             'in the genes argument. Possible values are: '
                             'hgnc_symbol, hgnc_id, and mgi_id.',
                        choices=['hgnc_symbol', 'hgnc_id', 'mgi_id'],
                        required=True)
    parser.add_argument('--stage', default='all',
                        help='The stage of processing to run.',
                        choices=['all', 'node_vectors', 'null_distribution',
                                 'statistics'])
    parser.add_argument('--base_folder', default=default_base_folder,
                        help='The base folder used to store GeneWalk '
                             'temporary and result files for a given project.')
    parser.add_argument('--network_source', default='pc',
                        help='The source of the network to be used.'
                             'Possible values are: pc, indra, edge_list, and '
                             'sif. In case of indra, edge_list, and sif, '
                             'the network_file argument must be specified.',
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
                             'multiprocessing environment.')
    parser.add_argument('--nreps', default=10, type=int,
                        help='The number of repeats to run when calculating '
                             'node vectors, and the null distribution.')
    parser.add_argument('--alpha_fdr', default=1, type=float,
                        help='The false discovery rate to use when '
                             'calculating the final statistics.')
    parser.add_argument('--save_dw', default=False, type=bool,
                        help='If True, the full DeepWalk object for each '
                             'repeat is saved in the project folder. This can '
                             'be useful for debugging but the files are '
                             'typically very large.')
    args = parser.parse_args()

    # Now we run the relevant stage of processing
    project_folder = create_project_folder(args.base_folder, args.project)
    if args.stage in ('all', 'node_vectors'):
        genes = read_gene_list(args.genes, args.id_type)
        save_pickle(genes, project_folder, 'genewalk_genes')
        MG = load_network(args.network_source, args.network_file, genes)
        save_pickle(MG.graph, project_folder, 'genewalk_mg')
        for i in range(args.nreps):
            logger.info('%s/%s' % (i + 1, args.nreps))
            DW = run_walk(MG.graph)

            # Pickle the node vectors (embeddings) and DW object
            if args.save_dw:
                save_pickle(DW, project_folder, 'genewalk_dw_%d' % i)
            nv = copy.deepcopy(DW.model.wv)
            save_pickle(nv, project_folder, 'genewalk_dw_nv_%d' % i)

    if args.stage in ('all', 'null_distribution'):
        MG = load_pickle(project_folder, 'genewalk_mg')
        srs = []
        for i in args.nreps:
            logger.info('%s/%s' % (i + 1, args.nreps))
            RG = get_rand_graph(MG)
            DW = run_walk(RG)

            # Pickle the node vectors (embeddings) and DW object
            if args.save_dw:
                save_pickle(DW, project_folder, 'genewalk_dw_rand_%d' % i)
            nv = copy.deepcopy(DW.model.wv)
            save_pickle(nv, project_folder, 'genewalk_dw_nv_rand_%d' % i)

            sr = get_null_distributions(RG, nv)
            srs.append(sr)
        srd = get_srd(srs)
        save_pickle(srd, project_folder, 'genewalk_dw_rand_simdists')

    if args.stage in ('all', 'statistics'):
        MG = load_pickle(project_folder, 'genewalk_mg')
        genes = load_pickle(project_folder, 'genewalk_genes')
        nvs = [load_pickle(project_folder, 'genewalk_dw_nv_rand_%d' % (i + 1))
               for i in range(args.nreps)]
        null_dist = load_pickle(project_folder, 'genewalk_dw_rand_simdists')
        GW = GeneWalk(MG, genes, nvs, null_dist)
        df = GW.generate_output(alpha_FDR=args.alpha_fdr)
        fname = os.path.join(project_folder, 'genewalk_results.csv')
        df.to_csv(fname, index=False)
