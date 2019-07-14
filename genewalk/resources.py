import os
import gzip
import shutil
import logging
import urllib.request

logger = logging.getLogger('genewalk.resources')


class ResourceManager(object):
    def __init__(self, base_folder=None):
        self.base_folder = base_folder if base_folder else \
            os.path.join(os.path.expanduser('~'), 'genewalk')
        self.resource_folder = self._get_resource_folder()
        logger.info('Using %s as resource folder.' % self.resource_folder)

    def get_go_obo(self):
        fname = os.path.join(self.resource_folder, 'go.obo')
        if not os.path.exists(fname):
            download_go(fname)
        return fname

    def get_goa_gaf(self):
        fname = os.path.join(self.resource_folder, 'goa_human.gaf')
        if not os.path.exists(fname):
            url_goa = ('http://geneontology.org/gene-associations/'
                       'goa_human.gaf.gz')
            download_gz(fname, url_goa, self.resource_folder)
        return fname

    def get_pc(self):
        fname = os.path.join(self.resource_folder,
                             'PathwayCommons11.All.hgnc.sif')
        if not os.path.exists(fname):
            url_pc = ('http://www.pathwaycommons.org/archives/PC2/v11/'
                      'PathwayCommons11.All.hgnc.sif.gz')
            download_gz(fname, url_pc, self.resource_folder)
        return fname

    def _get_resource_folder(self):
        resource_dir = os.path.join(self.base_folder, 'resources')

        if not os.path.isdir(resource_dir):
            try:
                os.makedirs(resource_dir)
            except Exception:
                logger.warning(resource_dir + ' already exists')
        return resource_dir

    def download_all(self):
        self.get_go_obo()
        self.get_goa_gaf()
        self.get_pc()


def download_go(fname):
    url = 'http://snapshot.geneontology.org/ontology/go.obo'
    logger.info('Downloading %s into %s' % (url, fname))
    urllib.request.urlretrieve(url, fname)


def download_gz(fname, url, resource_folder):
    logger.info('Downloading %s and extracting into %s' % (url, fname))
    gz_file = os.path.join(fname + '.gz')
    urllib.request.urlretrieve(url, gz_file)
    with gzip.open(gz_file, 'rb') as fin:
        with open(fname, 'wb') as fout:
            shutil.copyfileobj(fin, fout)


if __name__ == '__main__':
    # Download all the resources if this script is run directly
    ResourceManager().download_all()
