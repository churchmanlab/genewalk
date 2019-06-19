import os
import gzip
import shutil
import logging
import urllib.request
from . import __version__

logger = logging.getLogger(__name__)

home_dir = os.path.expanduser('~')
resource_dir = os.path.join(home_dir, '.genewalk', __version__)

if not os.path.isdir(resource_dir):
    try:
        os.makedirs(resource_dir)
    except Exception:
        logger.warning(resource_dir + ' already exists')


def download_go(fname):
    url = 'http://snapshot.geneontology.org/ontology/go.obo'
    logger.info('Downloading %s into %s' % (url, fname))
    urllib.request.urlretrieve(url, fname)


def download_goa(fname):
    url = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
    logger.info('Downloading %s and extracting into %s' % (url, fname))
    gz_file = os.path.join(resource_dir, 'goa_human.gaf.gz')
    urllib.request.urlretrieve(url, gz_file)
    with gzip.open(gz_file, 'rb') as fin:
        with open(fname, 'wb') as fout:
            shutil.copyfileobj(fin, fout)


def get_go_obo():
    fname = os.path.join(resource_dir, 'go.obo')
    if not os.path.exists(fname):
        download_go(fname)
    return fname


def get_goa_gaf():
    fname = os.path.join(resource_dir, 'goa_human.gaf')
    if not os.path.exists(fname):
        download_goa(fname)
    return fname
