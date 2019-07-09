import logging
from indra.databases import hgnc_client


logger = logging.getLogger('genewalk.gene_lists')


def _read_lines(fname):
    with open(fname, 'r') as fh:
        lines = [line.strip() for line in fh.readlines()]
    return lines


def load_hgnc_symbols(fname):
    lines = _read_lines(fname)
    refs = []
    for hgnc_symbol in lines:
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


def load_hgnc_ids(fname):
    lines = _read_lines(fname)
    refs = []
    for hgnc_id in lines:
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


def load_mgi_ids(fname):
    lines = _read_lines(fname)
    refs = []
    for mgi_id in lines:
        if mgi_id.startswith('MGI:'):
            mgi_id = mgi_id[4:]
        ref = {'HGNC_SYMBOL': None, 'HGNC': hgnc_id, 'UP': None,
               'MGI': mgi_id}
        hgnc_id = hgnc_client.get_hgnc_from_mouse(mgi_id)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for MGI ID %s' %
                           mgi_id)
            continue
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


