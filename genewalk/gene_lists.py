import re
import csv
import logging


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
    gene_mapper = GeneMapper(resource_manager)
    with open(fname, 'r') as fh:
        # This is to make the list unique while preserving
        # the original order as much as possible
        unique_lines = []
        for line in fh.readlines():
            line = line.strip()
            if line not in unique_lines:
                unique_lines.append(line)
    if id_type == 'hgnc_symbol':
        refs = map_hgnc_symbols(unique_lines, gene_mapper)
    elif id_type == 'hgnc_id':
        refs = map_hgnc_ids(unique_lines, gene_mapper)
    elif id_type == 'ensembl_id':
        refs = map_ensembl_ids(unique_lines, gene_mapper)
    elif id_type == 'mgi_id':
        refs = map_mgi_ids(unique_lines, gene_mapper)
    elif id_type == 'rgd_id':
        refs = map_rgd_ids(unique_lines, gene_mapper)
    elif id_type == 'entrez_human':
        refs = map_entrez_human(unique_lines, gene_mapper)
    elif id_type == 'entrez_mouse':
        refs = map_entrez_mouse(unique_lines, gene_mapper)
    elif id_type == 'custom':
        refs = [{'ID': c} for c in unique_lines]
    else:
        raise ValueError('Unknown id_type: %s' % id_type)
    if not refs:
        raise ValueError('None of the IDs in %s could be mapped. It is '
                         'likely that the file uses an ID type or format '
                         'that GeneWalk cannot interpret.' % fname)
    return refs


def map_hgnc_symbols(hgnc_symbols, gene_mapper):
    """Return references based on a list of HGNC symbols."""
    refs = []
    for hgnc_symbol in hgnc_symbols:
        ref = {'HGNC_SYMBOL': hgnc_symbol, 'HGNC': None, 'UP': None}
        hgnc_id = gene_mapper.get_current_hgnc_id(hgnc_symbol)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for symbol %s' % hgnc_symbol)
            continue
        elif isinstance(hgnc_id, list):
            logger.warning('More than one current HGNC ID for outdated '
                           'symbol %s' % hgnc_symbol)
            continue
        ref['HGNC'] = hgnc_id
        uniprot_id = gene_mapper.get_uniprot_id(hgnc_id)
        if not uniprot_id:
            logger.warning('Could not get UniProt ID for symbol %s' %
                           hgnc_symbol)
            continue
        ref['UP'] = uniprot_id
        refs.append(ref)
    return refs


def map_hgnc_ids(hgnc_ids, gene_mapper):
    """Return references based on a list of HGNC IDs."""
    refs = []
    for hgnc_id in hgnc_ids:
        if hgnc_id.startswith('HGNC:'):
            hgnc_id = hgnc_id[5:]
        hgnc_ref = _refs_from_hgnc_id(hgnc_id, gene_mapper)
        if hgnc_ref is None:
            continue
        refs.append(hgnc_ref)
    return refs


def _refs_from_hgnc_id(hgnc_id, gene_mapper):
    ref = {'HGNC_SYMBOL': None, 'HGNC': hgnc_id, 'UP': None}
    hgnc_name = gene_mapper.get_hgnc_name(hgnc_id)
    if not hgnc_name:
        logger.warning('Could not get HGNC name for ID %s' %
                       hgnc_id)
        return None
    ref['HGNC_SYMBOL'] = hgnc_name
    uniprot_id = gene_mapper.get_uniprot_id(hgnc_id)
    if not uniprot_id:
        logger.warning('Could not get UniProt ID for HGNC ID %s' %
                       hgnc_id)
        return None
    ref['UP'] = uniprot_id
    return ref


def map_mgi_ids(mgi_ids, gene_mapper):
    """Return references based on a list of MGI IDs."""
    refs = []
    for mgi_id in mgi_ids:
        if mgi_id.startswith('MGI:'):
            mgi_id = mgi_id[4:]
        mgi_ref = _refs_from_mgi_id(mgi_id, gene_mapper)
        if mgi_ref is None:
            continue
        refs.append(mgi_ref)
    return refs


def _refs_from_mgi_id(mgi_id, gene_mapper):
    ref = {'MGI': mgi_id}
    hgnc_id = gene_mapper.get_hgnc_from_mgi(mgi_id)
    if hgnc_id is None:
        logger.warning('Could not get HGNC ID for MGI ID %s' %
                       mgi_id)
        return None
    hgnc_ref = _refs_from_hgnc_id(hgnc_id, gene_mapper)
    if hgnc_ref is None:
        return None
    ref.update(hgnc_ref)
    return ref


def map_rgd_ids(rgd_ids, gene_mapper):
    """Return references based on a list of RGD IDs."""
    refs = []
    for rgd_id in rgd_ids:
        if rgd_id.startswith('RGD:'):
            rgd_id = rgd_id[4:]
        rgd_ref = _refs_from_rgd_id(rgd_id, gene_mapper)
        if rgd_ref is None:
            continue
        refs.append(rgd_ref)
    return refs


def _refs_from_rgd_id(rgd_id, gene_mapper):
    ref = {'RGD': rgd_id}
    hgnc_id = gene_mapper.get_hgnc_from_rgd(rgd_id)
    if hgnc_id is None:
        logger.warning('Could not get HGNC ID for RGD ID %s' %
                       rgd_id)
        return None
    hgnc_ref = _refs_from_hgnc_id(hgnc_id, gene_mapper)
    if hgnc_ref is None:
        return None
    ref.update(hgnc_ref)
    return ref


def map_ensembl_ids(ensembl_ids, gene_mapper):
    """Return references based on a list of Ensembl IDs."""
    refs = []
    for ensembl_id in ensembl_ids:
        ref = {'ENSEMBL': ensembl_id}
        ensembl_id = ensembl_id.split('.', maxsplit=1)[0]
        hgnc_id = gene_mapper.get_hgnc_from_ensembl(ensembl_id)
        if not hgnc_id:
            logger.warning('Could not get HGNC ID for ENSEMBL ID %s' %
                           ensembl_id)
            continue
        hgnc_ref = _refs_from_hgnc_id(hgnc_id, gene_mapper)
        if hgnc_ref is None:
            continue
        ref.update(hgnc_ref)
        refs.append(ref)
    return refs


def map_entrez_human(entrez_ids, gene_mapper):
    """Return references based on human Entrez gene IDs."""
    refs = []
    for entrez_id in entrez_ids:
        ref = {'EGID': entrez_id}
        hgnc_id = gene_mapper.get_hgnc_from_entrez(entrez_id)
        if hgnc_id is None:
            logger.warning("Could not find HGNC ID for Entrez ID %s" %
                           entrez_id)
            continue
        hgnc_ref = _refs_from_hgnc_id(hgnc_id, gene_mapper)
        if hgnc_ref is None:
            continue
        ref.update(hgnc_ref)
        refs.append(ref)
    return refs


def map_entrez_mouse(entrez_ids, gene_mapper):
    """Return references based on mouse Entrez gene IDs."""
    # Get the entrez file path from the resource manager
    refs = []
    for entrez_id in entrez_ids:
        mgi_id = gene_mapper.entrez_to_mgi.get(entrez_id)
        if not mgi_id:
            logger.warning("Could not find an MGI mapping for Entrez ID %s"
                           % entrez_id)
            continue
        ref = {'EGID': entrez_id, 'MGI': mgi_id}
        mgi_refs = _refs_from_mgi_id(mgi_id, gene_mapper)
        if mgi_refs is None:
            continue
        ref.update(mgi_refs)
        refs.append(ref)
    return refs


class GeneMapper:
    def __init__(self, resource_manager):
        self.resource_manager = resource_manager
        self.hgnc_file = self.resource_manager.get_hgnc()
        self.mgi_entrez_file = self.resource_manager.get_mgi_entrez()

        # Process the MGI-Entrez mapping file
        self.entrez_to_mgi = {}
        with open(self.mgi_entrez_file, 'r') as fh:
            csvreader = csv.reader(fh, delimiter='\t')
            for row in csvreader:
                # Remove "MGI:" prefix
                mgi = row[0][4:]
                entrez = row[8]
                self.entrez_to_mgi[entrez] = mgi

        self.hgnc_id_to_name = {}
        self.hgnc_name_to_id = {}
        self.hgnc_withdrawn_to_new = {}
        self.hgnc_to_uniprot = {}
        self.mgi_to_hgnc = {}
        self.rgd_to_hgnc = {}
        self.entrez_to_hgnc = {}
        self.ensembl_to_hgnc = {}
        self.prev_sym_map = {}

        with open(self.hgnc_file, 'r', encoding='utf-8') as fh:
            csvreader = csv.reader(fh, delimiter='\t')
            # Skip the header
            next(csvreader)
            for row in csvreader:
                hgnc_id, hgnc_name, description, prev_sym_entry, hgnc_status,\
                    entrez_id, uniprot_id, mgi_id, rgd_id, ensembl_id = row
                hgnc_id = hgnc_id[5:]
                if hgnc_status in {'Approved', 'Entry Withdrawn'}:
                    self.hgnc_id_to_name[hgnc_id] = hgnc_name
                    # Note that withdrawn entries don't overlap with approved
                    # entries at this point so it's safe to add mappings for
                    # withdrawn names
                    self.hgnc_name_to_id[hgnc_name] = hgnc_id
                elif hgnc_status == 'Symbol Withdrawn':
                    m = re.match(r'symbol withdrawn, see \[HGNC:(?: ?)(\d+)\]',
                                 description)
                    new_id = m.groups()[0]
                    self.hgnc_withdrawn_to_new[hgnc_id] = new_id
                # Uniprot
                if uniprot_id:
                    self.hgnc_to_uniprot[hgnc_id] = uniprot_id
                # Entrez
                if entrez_id:
                    self.entrez_to_hgnc[entrez_id] = hgnc_id
                # Mouse
                if mgi_id:
                    mgi_ids = mgi_id.split(', ')
                    for mgi_id in mgi_ids:
                        if mgi_id.startswith('MGI:'):
                            mgi_id = mgi_id[4:]
                        self.mgi_to_hgnc[mgi_id] = hgnc_id
                # Rat
                if rgd_id:
                    rgd_ids = rgd_id.split(', ')
                    for rgd_id in rgd_ids:
                        if rgd_id.startswith('RGD:'):
                            rgd_id = rgd_id[4:]
                        self.rgd_to_hgnc[rgd_id] = hgnc_id
                # Previous symbols
                if prev_sym_entry:
                    prev_syms = prev_sym_entry.split(', ')
                    for prev_sym in prev_syms:
                        # If we already mapped this previous symbol
                        # to another ID
                        if prev_sym in self.prev_sym_map:
                            # If we already have a list here, we just extend it
                            if isinstance(self.prev_sym_map[prev_sym], list):
                                self.prev_sym_map[prev_sym].append(hgnc_id)
                            # Otherwise we create a list and start it with the
                            # two IDs we know the symbol is mapped to
                            else:
                                self.prev_sym_map[prev_sym] = \
                                    [self.prev_sym_map[prev_sym], hgnc_id]
                        # Otherwise we just make a string entry here
                        else:
                            self.prev_sym_map[prev_sym] = hgnc_id
                # Ensembl IDs
                if ensembl_id:
                    self.ensembl_to_hgnc[ensembl_id] = hgnc_id
            for old_id, new_id in self.hgnc_withdrawn_to_new.items():
                self.hgnc_id_to_name[old_id] = self.hgnc_id_to_name[new_id]

    def get_hgnc_name(self, hgnc_id):
        """Return the HGNC symbol corresponding to the given HGNC ID.

        Parameters
        ----------
        hgnc_id : str
            The HGNC ID to be converted.

        Returns
        -------
        hgnc_name : str
            The HGNC symbol corresponding to the given HGNC ID.
        """
        hgnc_name = self.hgnc_id_to_name.get(hgnc_id)
        return hgnc_name

    def get_hgnc_id(self, hgnc_name):
        """Return the HGNC ID corresponding to the given HGNC symbol.

        Parameters
        ----------
        hgnc_name : str
            The HGNC symbol to be converted. Example: BRAF

        Returns
        -------
        hgnc_id : str
            The HGNC ID corresponding to the given HGNC symbol.
        """
        return self.hgnc_name_to_id.get(hgnc_name)

    def get_current_hgnc_id(self, hgnc_name):
        """Return HGNC ID(s) corresponding to a current or outdated HGNC symbol.

        Parameters
        ----------
        hgnc_name : str
            The HGNC symbol to be converted, possibly an outdated symbol.

        Returns
        -------
        str or list of str or None
            If there is a single HGNC ID corresponding to the given current or
            outdated HGNC symbol, that ID is returned as a string. If the symbol
            is outdated and maps to multiple current IDs, a list of these
            IDs is returned. If the given name doesn't correspond to either
            a current or an outdated HGNC symbol, None is returned.
        """
        hgnc_id = self.get_hgnc_id(hgnc_name)
        if hgnc_id:
            return hgnc_id
        hgnc_id = self.prev_sym_map.get(hgnc_name)
        return hgnc_id

    def get_uniprot_id(self, hgnc_id):
        """Return the UniProt ID corresponding to the given HGNC ID.

        Parameters
        ----------
        hgnc_id : str
            The HGNC ID to be converted. Note that the HGNC ID is a number that is
            passed as a string. It is not the same as the HGNC gene symbol.

        Returns
        -------
        uniprot_id : str
            The UniProt ID corresponding to the given HGNC ID.
        """
        uniprot_id = self.hgnc_to_uniprot.get(hgnc_id)
        # The lookup can yield an empty string. Instead return None.
        if not uniprot_id:
            return None
        return uniprot_id

    def get_hgnc_from_entrez(self, entrez_id):
        """Return the HGNC ID corresponding to the given Entrez ID.

        Parameters
        ----------
        entrez_id : str
            The Entrez ID to be converted, a number passed as a string.

        Returns
        -------
        hgnc_id : str
            The HGNC ID corresponding to the given Entrez ID.
        """
        hgnc_id = self.entrez_to_hgnc.get(entrez_id)
        return hgnc_id

    def get_hgnc_from_ensembl(self, ensembl_id):
        """Return the HGNC ID corresponding to the given Ensembl ID.

        Parameters
        ----------
        ensembl_id : str
            The Ensembl ID to be converted, a number passed as a string.

        Returns
        -------
        hgnc_id : str
            The HGNC ID corresponding to the given Ensembl ID.
        """
        return self.ensembl_to_hgnc.get(ensembl_id)

    def get_hgnc_from_mgi(self, mgi_id):
        """Return the HGNC ID corresponding to the given MGI mouse gene ID.

        Parameters
        ----------
        mgi_id : str
            The MGI ID to be converted. Example: "2444934"

        Returns
        -------
        hgnc_id : str
            The HGNC ID corresponding to the given MGI ID.
        """
        if mgi_id.startswith('MGI:'):
            mgi_id = mgi_id[4:]
        return self.mgi_to_hgnc.get(mgi_id)

    def get_hgnc_from_rgd(self, rgd_id):
        """Return the HGNC ID corresponding to the given RGD rat gene ID.

        Parameters
        ----------
        rgd_id : str
            The RGD ID to be converted. Example: "1564928"

        Returns
        -------
        hgnc_id : str
            The HGNC ID corresponding to the given RGD ID.
        """
        if rgd_id.startswith('RGD:'):
            rgd_id = rgd_id[4:]
        return self.rgd_to_hgnc.get(rgd_id)
