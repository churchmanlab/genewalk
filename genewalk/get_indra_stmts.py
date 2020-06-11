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
from indra.util import batch_iter
from indra.databases import go_client
from indra.sources import indra_db_rest
from indra.databases import hgnc_client
from indra.ontology.bio import bio_ontology


logger = logging.getLogger('genewalk.get_indra_stmts')


def load_genes(fname):
    """Return a list of genes IDs from a file with lines like HGNC:123."""
    with open(fname, 'r') as fh:
        # Get the HGNC IDs from the list of genes file, assuming that
        # each line looks like HGNC:123
        genes = [l.strip().split(':')[1] for l in fh.readlines()]
    logger.info('Loaded %d genes from %s' % (len(genes), fname))
    return genes


def load_mouse_genes(fname):
    """Return a list of human genes based on a table of mouse genes."""
    # assumes the csv has headers
    df = pandas.read_csv(fname)
    for c in df.columns:
        # assumes the first column starting with MGI is the relevant one
        # with MGI:IDs
        if c.startswith('MGI'):
            df = df.rename(columns={c: 'MGI'})
            break
    mgi_ids = df['MGI']
    genes = []
    for mgi_id in mgi_ids:
        if mgi_id.startswith('MGI:'):
            mgi_id = mgi_id[4:]
        hgnc_id = hgnc_client.get_hgnc_from_mouse(mgi_id)
        if not hgnc_id:
            print('Could not find human gene corresponding to MGI %s' % mgi_id)
            continue
        genes.append(hgnc_id)
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


def get_famplex_terms(genes):
    """Get a list of associated FamPlex IDs from a list of gene IDs."""
    all_parents = set()
    for hgnc_id in genes:
        parent_ids = {p[1] for p in
                      bio_ontology.get_parents('HGNC', hgnc_id)}
        all_parents |= parent_ids
    fplx_terms = sorted(list(all_parents))
    logger.info('Found %d relevant FamPlex terms.' % (len(fplx_terms)))
    return fplx_terms


def get_famplex_links(df, fname):
    """Given a list of INDRA Statements, construct FamPlex links."""
    genes_appearing = (set(df[df.agA_ns == 'HGNC'].agA_name) |
                       set(df[df.agB_ns == 'HGNC'].agB_name))
    fplx_appearing = (set(df[df.agA_ns == 'FPLX'].agA_id) |
                      set(df[df.agB_ns == 'FPLX'].agB_id))
    links = get_famplex_links_from_lists(genes_appearing, fplx_appearing)
    with open(fname, 'w') as fh:
        for link in links:
            fh.write('%s,%s\n' % link)


def get_famplex_links_from_stmts(stmts):
    genes_appearing = set()
    fplx_appearing = set()
    for stmt in stmts:
        agents = [a for a in stmt.agent_list() if a is not None]
        if len(agents) < 2:
            continue
        for agent in agents:
            if 'HGNC' in agent.db_refs:
                genes_appearing.add(agent.name)
            elif 'FPLX' in agent.db_refs:
                fplx_appearing.add(agent.name)
    return get_famplex_links_from_lists(genes_appearing, fplx_appearing)


def get_famplex_links_from_lists(genes_appearing, fplx_appearing):
    links = []
    for gene in genes_appearing:
        parent_ids = [p[1] for p in bio_ontology.get_parents('HGNC', gene)]
        parents_appearing = fplx_appearing & set(parent_ids)
        links += [(gene, parent) for parent in parents_appearing]
    for fplx_child in fplx_appearing:
        parent_ids = [p[1] for p in
                      bio_ontology.get_parents('FPLX', fplx_child)]
        parents_appearing = fplx_appearing & set(parent_ids)
        links += [(fplx_child, parent) for parent in parents_appearing]
    return links


def download_statements(df, ev_limit=5):
    """Download the INDRA Statements corresponding to entries in a data frame.
    """
    all_stmts = []
    for idx, group in enumerate(batch_iter(df.hash, 500)):
        logger.info('Getting statement batch %d' % idx)
        idbp = indra_db_rest.get_statements_by_hash(list(group),
                                                    ev_limit=ev_limit)
        all_stmts += idbp.statements
    return all_stmts


def remap_go_ids(stmts):
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None and 'GO' in agent.db_refs:
                prim_id = go_client.get_primary_id(agent.db_refs['GO'])
                if prim_id:
                    agent.db_refs['GO'] = prim_id


if __name__ == '__main__':
    # Handle command line arguments
    # TODO: make specific sets of statements for the use cases available in
    #  the repo (or on S3) and parameterize here to be able to load them.
    parser = argparse.ArgumentParser(
        description='Choose a file with a list of genes to get a SIF for.')
    parser.add_argument('--df', default='data/stmt_df.pkl')
    parser.add_argument('--genes', default='data/JQ1_HGNCidForINDRA.csv')
    parser.add_argument('--mouse_genes')
    parser.add_argument('--stmts', default='data/JQ1_HGNCidForINDRA_stmts.pkl')
    args = parser.parse_args()

    # Load genes and get FamPlex terms
    if args.mouse_genes:
        genes = load_mouse_genes(args.mouse_genes)
    else:
        genes = load_genes(args.genes)
    fplx_terms = get_famplex_terms(genes)
    # Load INDRA Statements in a flat data frame
    df = load_indra_df(args.df)
    # Filter the data frame to relevant entities
    df = filter_to_genes(df, genes, fplx_terms)
    # Download the Statement corresponding to each row
    stmts = download_statements(df)
    # Remap any outdated GO IDs
    remap_go_ids(stmts)
    # Dump the Statements into a pickle file
    dump_pickle(stmts, args.stmts)
