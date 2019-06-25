import os
import gzip
import shutil
import logging
import urllib.request

logger = logging.getLogger('genewalk.resources')

home_dir = os.path.expanduser('/n/groups/churchman/ri23/')#TEMP###################, change back to '~'
resource_dir = os.path.join(home_dir, 'genewalk','resources')

logger.addHandler(logging.FileHandler(os.path.join(resource_dir,'%s.log' % logger.name))) 

if not os.path.isdir(resource_dir):
    try:
        os.makedirs(resource_dir)
    except Exception:
        logger.warning(resource_dir + ' already exists')


def download_go(fname):
    url = 'http://snapshot.geneontology.org/ontology/go.obo'
    logger.info('Downloading %s into %s' % (url, fname))
    urllib.request.urlretrieve(url, fname)


def download_gz(fname,url):
    logger.info('Downloading %s and extracting into %s' % (url, fname))
    gz_file = os.path.join(resource_dir, fname+'.gz')
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
        url_goa = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
        download_gz(fname,url_goa)
    return fname


def get_pc():
    fname = os.path.join(resource_dir, 'PathwayCommons11.All.hgnc.sif')
    if not os.path.exists(fname):
        url_pc = 'http://www.pathwaycommons.org/archives/PC2/v11/PathwayCommons11.All.hgnc.sif.gz'
        download_gz(fname,url_pc)
    return fname