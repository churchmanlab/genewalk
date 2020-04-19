import csv
import pickle
import logging
from io import StringIO
from os.path import join, abspath, dirname
from collections import Counter
import requests
from indra.databases import hgnc_client, uniprot_client


logger = logging.getLogger('genewalk.gene_lists')


def read_gene_list(fname, id_type, resource_manager):
    """Return references for genes from a file with the given ID type.

    Parameters
    ----------
    fname : str
        The name of the file containing the list of genes. Each line of the
        file corresponds to a single gene.
    id_type : str
        The type of identifier contained in each line of the gene list file.
        Possible values are: hgnc_symbol, hgnc_id, ensembl_id, mgi_id.
    resource_manager : genewalk.resources.ResourceManager
        ResourceManager object, used to obtain entrez-mgi mappings if
        necessary.

    Returns
    -------
    dict
        A dictionary of references with keys including HGNCSYMBOL, HGNC, UP,
        and if id_type is mgi_id, MGI, with values corresponding to the
        identifiers of the provided list of genes.
    """
    with open(fname, 'r') as fh:
        # This is to make the list unique while preserving
        # the original order as much as possible
        unique_lines = []
        for line in fh.readlines():
            line = line.strip()
            if line not in unique_lines:
                unique_lines.append(line)
    if id_type == 'hgnc_symbol':
        refs = map_hgnc_symbols(unique_lines)
    elif id_type == 'hgnc_id':
        refs = map_hgnc_ids(unique_lines)
    elif id_type == 'ensembl_id':
        refs = map_ensembl_ids(unique_lines)
    elif id_type == 'mgi_id':
        refs = map_mgi_ids(unique_lines)
    elif id_type == 'entrez_human':
        refs = map_entrez_human(unique_lines)
    elif id_type == 'entrez_mouse':
        refs = map_entrez_mouse(unique_lines, resource_manager)
    else:
        raise ValueError('Unknown id_type: %s' % id_type)
    if not refs:
        raise ValueError('None of the IDs in %s could be mapped. It is '
                         'likely that the file uses an ID type or format '
                         'that GeneWalk cannot interpret.' % fname)
    return refs


def map_hgnc_symbols(hgnc_symbols):
    """Return references based on a list of HGNC symbols."""
    refs = []
    for hgnc_symbol in hgnc_symbols:
        ref = {'HGNC_SYMBOL': hgnc_symbol, 'HGNC': None, 'UP': None}
        hgnc_id = hgnc_client.get_current_hgnc_id(hgnc_symbol)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for symbol %s' % hgnc_symbol)
            continue
        elif isinstance(hgnc_id, list):
            logger.warning('More than one current HGNC ID for outdated '
                           'symbol %s' % hgnc_symbol)
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
        hgnc_ref = _refs_from_hgnc_id(hgnc_id)
        if hgnc_ref is None:
            continue
        refs.append(hgnc_ref)
    return refs


def _refs_from_hgnc_id(hgnc_id):
    ref = {'HGNC_SYMBOL': None, 'HGNC': hgnc_id, 'UP': None}
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    if not hgnc_name:
        logger.warning('Could not get HGNC name for ID %s' %
                       hgnc_id)
        return None
    ref['HGNC_SYMBOL'] = hgnc_name
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    if not uniprot_id:
        logger.warning('Could not get UniProt ID for HGNC ID %s' %
                       hgnc_id)
        return None
    ref['UP'] = uniprot_id
    return ref


def map_mgi_ids(mgi_ids):
    """Return references based on a list of MGI IDs."""
    refs = []
    for mgi_id in mgi_ids:
        if mgi_id.startswith('MGI:'):
            mgi_id = mgi_id[4:]
        mgi_ref = _refs_from_mgi_id(mgi_id)
        if mgi_ref is None:
            continue
        refs.append(mgi_ref)
    return refs


def _refs_from_mgi_id(mgi_id):
    ref = {'MGI': mgi_id}
    hgnc_id = hgnc_client.get_hgnc_from_mouse(mgi_id)
    if hgnc_id is None:
        logger.warning('Could not get HGNC ID for MGI ID %s' %
                       mgi_id)
        return None
    hgnc_ref = _refs_from_hgnc_id(hgnc_id)
    if hgnc_ref is None:
        return None
    ref.update(hgnc_ref)
    return ref


def map_ensembl_ids(ensembl_ids):
    """Return references based on a list of Ensembl IDs."""
    refs = []
    for ensembl_id in ensembl_ids:
        ref = {'ENSEMBL': ensembl_id}
        ensembl_id = ensembl_id.split('.', maxsplit=1)[0]
        hgnc_id = hgnc_client.get_hgnc_from_ensembl(ensembl_id)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for ENSEMBL ID %s' %
                           ensembl_id)
            continue
        hgnc_ref = _refs_from_hgnc_id(hgnc_id)
        if hgnc_ref is None:
            continue
        ref.update(hgnc_ref)
        refs.append(ref)
    return refs


def map_entrez_human(entrez_ids):
    """Return references based on human Entrez gene IDs."""
    refs = []
    for entrez_id in entrez_ids:
        ref = {'EGID': entrez_id}
        hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
        if hgnc_id is None:
            logger.warning("Could not find HGNC ID for Entrez ID %s" %
                           entrez_id)
            continue
        hgnc_ref = _refs_from_hgnc_id(hgnc_id)
        if hgnc_ref is None:
            continue
        ref.update(hgnc_ref)
        refs.append(ref)
    return refs


def map_entrez_mouse(entrez_ids, rm):
    """Return references based on mouse Entrez gene IDs."""
    # Get the entrez file path from the resource manager
    mgi_entrez_file = rm.get_mgi_entrez()
    # Process the MGI-Entrez mapping file
    entrez_to_mgi = {}
    with open(mgi_entrez_file, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for row in csvreader:
            # Remove "MGI:" prefix
            mgi = row[0][4:]
            entrez = row[8]
            entrez_to_mgi[entrez] = mgi
    refs = []
    for entrez_id in entrez_ids:
        mgi_id = entrez_to_mgi.get(entrez_id)
        if not mgi_id:
            logger.warning("Could not find an MGI mapping for Entrez ID %s"
                           % entrez_id)
            continue
        ref = {'EGID': entrez_id, 'MGI': mgi_id}
        mgi_refs = _refs_from_mgi_id(mgi_id)
        if mgi_refs is None:
            continue
        ref.update(mgi_refs)
        refs.append(ref)
    return refs
