import os
import glob
import shutil
import logging
from genewalk.cli import default_base_folder

logger = logging.getLogger(__name__)

TEST_MODULE = os.path.dirname(os.path.abspath(__file__))
TEST_RESOURCES = os.path.join(TEST_MODULE, 'resources')
TEST_BASE_FOLDER = os.path.join(default_base_folder, '.test')


def place_resource_files():
    test_resource_files = glob.glob(os.path.join(TEST_RESOURCES, '*'))
    test_resource_folder = \
        os.path.join(TEST_BASE_FOLDER, 'resources')
    os.makedirs(test_resource_folder, exist_ok=True)
    for test_file in test_resource_files:
        logger.debug('Copying %s into %s' % (test_file, test_resource_folder))
        shutil.copy(test_file,
                    os.path.join(test_resource_folder,
                                 os.path.basename(test_file)))
