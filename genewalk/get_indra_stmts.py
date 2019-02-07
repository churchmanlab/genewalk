"""This script creates a pickle file of INDRA Statements for a specific set of
genes, based on a pre-assembled network of interactions produced by INDRA. It
loads the INDRA interactions as a pandas DataFrame and, filters it to the genes
and familiex/complexes of interest, as well as targeted biological processes.
It downloads the relevant Statement objects for reference and dumps them into
a pickle file.
"""
import pandas
import pickle
import logging
import argparse
import itertools
from indra.util import batch_iter
from indra.sources import indra_db_rest
from indra.databases import hgnc_client
from indra.preassembler.hierarchy_manager import hierarchies


logger = logging.getLogger('genewalk.make_indra_sif')


def load_genes(fname):
    """Return a list of genes IDs from a file with lines like HGNC:123."""
    with open(fname, 'r') as fh:
        # Get the HGNC IDs from the list of genes file, assuming that
        # each line looks like HGNC:123
        genes = [l.strip().split(':')[1] for l in fh.readlines()]
    logger.info('Loaded %d genes from %s' % (len(genes), fname))
    return genes


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def dump_pickle(stmts, fname):
    """Dump a list of Statements into a picke file."""
    with open(fname, 'wb') as fh:
        pickle.dump(stmts, fh)
    logger.info('Dumped %d statements into %s' % (len(stmts), fname))


def filter_to_genes(df, genes, fplx_terms):
    """Filter a data frame of INDRA Statements given gene and FamPlex IDs."""
    # Look for sources that are in the gene list or whose families/complexes
    # are in the FamPlex term list
    source_filter = (((df.agA_ns == 'HGNC') & (df.agA_id.isin(genes))) |
                     ((df.agA_ns == 'FPLX') & (df.agA_id.isin(fplx_terms))))
    # Look for targets that are in the gene list or whose families/complexes
    # are in the FamPlex term list, or which are GO terms
    target_filter = (((df.agB_ns == 'HGNC') & (df.agB_id.isin(genes))) |
                     ((df.agB_ns == 'FPLX') & (df.agB_id.isin(fplx_terms))) |
                     (df.agB_ns == 'GO'))
    # sources/targets
    df = df[source_filter & target_filter]
    logger.info('Filtered data frame to %d rows.' % len(df))
    return df


def get_gene_parents(hgnc_name):
    eh = hierarchies['entity']
    gene_uri = eh.get_uri('HGNC', hgnc_name)
    parents = eh.get_parents(gene_uri)
    parent_ids = [eh.ns_id_from_uri(par_uri)[1] for par_uri in parents]
    return parent_ids


def get_famplex_terms(genes):
    """Get a list of associated FamPlex IDs from a list of gene IDs."""
    all_parents = set()
    for hgnc_id in genes:
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        parent_ids = get_gene_parents(hgnc_name)
        all_parents |= set(parent_ids)
    fplx_terms = sorted(list(all_parents))
    logger.info('Found %d relevant FamPlex terms.' % (len(fplx_terms)))
    return fplx_terms


def get_famplex_links(df, fname):
    """Given a list of INDRA Statements, construct FamPlex links."""
    genes_appearing = (set(df[df.agA_ns == 'HGNC'].agA_name) |
                       set(df[df.agB_ns == 'HGNC'].agB_name))
    fplx_appearing = (set(df[df.agA_ns == 'FPLX'].agA_id) |
                      set(df[df.agB_ns == 'FPLX'].agB_id))
    links = []
    for gene in genes_appearing:
        parent_ids = get_gene_parents(gene)
        parents_appearing = fplx_appearing & set(parent_ids)
        links += [(gene, parent) for parent in parents_appearing]
    eh = hierarchies['entity']
    for fplx_child in fplx_appearing:
        fplx_uri = eh.get_uri('FPLX', fplx_child)
        parents = eh.get_parents(fplx_uri)
        parent_ids = [eh.ns_id_from_uri(par_uri)[1] for par_uri in parents]
        parents_appearing = fplx_appearing & set(parent_ids)
        links += [(fplx_child, parent) for parent in parents_appearing]
    with open(fname, 'w') as fh:
        for link in links:
            fh.write('%s,%s\n' % link)
    return links


def download_statements(df):
    """Download the INDRA Statements corresponding to entries in a data frame.
    """
    from indra.sources.indra_db_rest.util import logger
    logger.setLevel(logging.ERROR)
    all_stmts = []
    for idx, group in enumerate(batch_iter(df.hash, 500)):
        logger.info('Getting statement batch %d' % idx)
        stmts = indra_db_rest.get_statements_by_hash(list(group))
        all_stmts += stmts
    return all_stmts


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description='Choose a file with a list of genes to get a SIF for.')
    parser.add_argument('--df', default='data/stmt_df.pkl')
    parser.add_argument('--genes', default='data/JQ1_HGNCidForINDRA.csv')
    parser.add_argument('--stmts', default='data/JQ1_HGNCidForINDRA.pkl')
    parser.add_argument('--fplx', default='data/JQ1_HGNCidForINDRA.csv')
    args = parser.parse_args()
    # Load genes and get FamPlex terms
    genes = load_genes(args.genes)
    fplx_terms = get_famplex_terms(genes)
    # Load INDRA Statements in a flat data frame
    df = load_indra_df(args.df)
    # Filter the data frame to relevant entities
    df = filter_to_genes(df, genes, fplx_terms)
    # Download the Statement corresponding to each row
    # stmts = download_statements(df)
    # Dump the Statements into a pickle file
    # dump_pickle(stmts, args.stmts)
    fplx_links = get_famplex_links(df, args.fplx)
