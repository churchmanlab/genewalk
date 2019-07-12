import logging
from indra.databases import hgnc_client


logger = logging.getLogger('genewalk.gene_lists')

# TODO: map to MGI symbols if the original genes were mouse genes


def read_gene_list(fname, id_type):
    """Return references for genes from a file with the given ID type.

    Parameters
    ----------
    fname : str
        The name of the file containing the list of genes. Each line of the
        file corresponds to a single gene.
    id_type : str
        The type of identifier contained in each line of the gene list file.
        Possible values are: hgnc_symbol, hgnc_id, mgi_id.

    Returns
    -------
    dict
        A dictionary of references with keys including HGNCSYMBOL, HGNC, UP,
        and if id_type is mgi_id, MGI, with values corresponding to the
        identifiers of the provided list of genes.
    """
    with open(fname, 'r') as fh:
        lines = [line.strip() for line in fh.readlines()]
    if id_type == 'hgnc_symbol':
        return map_hgnc_symbols(lines)
    elif id_type == 'hgnc_id':
        return map_hgnc_ids(lines)
    elif id_type == 'mgi_id':
        return map_mgi_ids(lines)
    else:
        raise ValueError('Unknown id_type: %s' % id_type)


def map_hgnc_symbols(hgnc_symbols):
    """Return references based on a list of HGNC symbols."""
    refs = []
    for hgnc_symbol in hgnc_symbols:
        ref = {'HGNC_SYMBOL': hgnc_symbol, 'HGNC': None, 'UP': None}
        hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for symbol %s' % hgnc_symbol)
            continue
        ref['HGNC'] = hgnc_id
        uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
        if not uniprot_id:
            logger.warning('Could not get UniProt ID for symbol %s' %
                           hgnc_symbol)
            continue
        ref['UP'] = uniprot_id
        refs.append(ref)
    return refs


def map_hgnc_ids(hgnc_ids):
    """Return references based on a list of HGNC IDs."""
    refs = []
    for hgnc_id in hgnc_ids:
        if hgnc_id.startswith('HGNC:'):
            hgnc_id = hgnc_id[5:]
        ref = {'HGNC_SYMBOL': None, 'HGNC': hgnc_id, 'UP': None}
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        if not hgnc_name:
            logger.warning('Could not get HGNC name for ID %s' %
                           hgnc_id)
            continue
        ref['HGNC_SYMBOL'] = hgnc_name
        uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
        if not uniprot_id:
            logger.warning('Could not get UniProt ID for HGNC ID %s' %
                           hgnc_id)
            continue
        ref['UP'] = uniprot_id
        refs.append(ref)
    return refs


def map_mgi_ids(mgi_ids):
    """Return references based on a list of MGI IDs."""
    refs = []
    for mgi_id in mgi_ids:
        if mgi_id.startswith('MGI:'):
            mgi_id = mgi_id[4:]
        ref = {'HGNC_SYMBOL': None, 'HGNC': None, 'UP': None,
               'MGI': mgi_id}
        hgnc_id = hgnc_client.get_hgnc_from_mouse(mgi_id)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for MGI ID %s' %
                           mgi_id)
            continue
        ref['HGNC'] = hgnc_id
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        if not hgnc_name:
            logger.warning('Could not get HGNC name for ID %s' %
                           hgnc_id)
            continue
        ref['HGNC_SYMBOL'] = hgnc_name
        uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
        if not uniprot_id:
            logger.warning('Could not get UniProt ID for HGNC ID %s' %
                           hgnc_id)
            continue
        ref['UP'] = uniprot_id
        refs.append(ref)
    return refs
