import os
import gzip
import shutil
import logging
import urllib.request
import pandas
from indra.databases import hgnc_client

logger = logging.getLogger('genewalk.resources')


class ResourceManager(object):
    """Class to manage the download, caching and access of resource files."""
    def __init__(self, base_folder=None):
        self.base_folder = base_folder if base_folder else \
            os.path.join(os.path.expanduser('~'), 'genewalk')
        self.resource_folder = self._get_resource_folder()
        logger.info('Using %s as resource folder.' % self.resource_folder)

    def get_go_obo(self):
        fname = os.path.join(self.resource_folder, 'go.obo')
        if not os.path.exists(fname):
            url = 'http://snapshot.geneontology.org/ontology/go.obo'
            download_url(url, fname)
        return fname

    def get_goa_gaf(self):
        fname = os.path.join(self.resource_folder, 'goa_human.gaf')
        if not os.path.exists(fname):
            url_goa = ('http://geneontology.org/gene-associations/'
                       'goa_human.gaf.gz')
            download_gz(fname, url_goa)
        return fname

    def get_pc(self):
        fname_current = os.path.join(self.resource_folder,
                                     'PathwayCommons12.All.hgnc_current.sif')
        if not os.path.exists(fname_current):
            fname = os.path.join(self.resource_folder,
                                 'PathwayCommons12.All.hgnc.sif')
            if not os.path.exists(fname):
                url_pc = ('http://www.pathwaycommons.org/archives/PC2/v12/'
                          'PathwayCommons12.All.hgnc.sif.gz')
                download_gz(fname, url_pc)
            self._replace_outdated_hgnc_symbols(fname,fname_current)
        return fname_current

    def get_mgi_entrez(self):
        fname = os.path.join(self.resource_folder, 'MGI_EntrezGene.rpt')
        if not os.path.exists(fname):
            url = 'http://www.informatics.jax.org/downloads/reports/' \
                  'MGI_EntrezGene.rpt'
            download_url(url, fname)
        return fname

    def _get_resource_folder(self):
        resource_dir = os.path.join(self.base_folder, 'resources')
        if not os.path.isdir(resource_dir):
            try:
                os.makedirs(resource_dir)
            except Exception:
                logger.warning(resource_dir + ' already exists')
        return resource_dir

    def _replace_outdated_hgnc_symbols(self,pc_old,pc_current):
        logger.info('Replacing outdated HGNC symbols in %s and save as %s' % \
                    (pc_old, pc_current))
        pc = pandas.read_csv(pc_old,sep='\t',dtype=str, header=None)
        col_mapper = {}
        col_mapper[0] = 'source'
        col_mapper[1] = 'rel_type'
        col_mapper[2] = 'target'
        pc = pc.rename(mapper=col_mapper, axis='columns')
        all_symbols = set(pc['source']).union(pc['target'])
        symbol_map = {}
        for sym in all_symbols:
            if not sym.startswith('CHEBI:'):
                hgnc_id = hgnc_client.get_current_hgnc_id(sym)
                if not hgnc_id:
                    continue
                elif isinstance(hgnc_id, list):
                    #outdated gene symbol is ambiguous: maps to multiple genes
                    continue
                latest_symbol = hgnc_client.get_hgnc_name(hgnc_id)
                if latest_symbol != sym:
                    symbol_map[sym] = latest_symbol
        if symbol_map:
            pc.replace(symbol_map,inplace=True)         
        pc.to_csv(pc_current, sep='\t', header=False, index=False)
        os.remove(pc_old)
        
    def download_all(self):
        self.get_go_obo()
        self.get_goa_gaf()
        self.get_pc()
        self.get_mgi_entrez()
        

def download_url(url, fname):
    logger.info('Downloading %s into %s' % (url, fname))
    urllib.request.urlretrieve(url, fname)


def download_gz(fname, url):
    logger.info('Downloading %s and extracting into %s' % (url, fname))
    gz_file = fname + '.gz'
    urllib.request.urlretrieve(url, gz_file)
    with gzip.open(gz_file, 'rb') as fin:
        with open(fname, 'wb') as fout:
            shutil.copyfileobj(fin, fout)


if __name__ == '__main__':
    # Download all the resources if this script is run directly
    ResourceManager().download_all()
