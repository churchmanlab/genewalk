import logging
default_logger_format = ('%(levelname)s: [%(asctime)s] %(name)s'
                         ' - %(message)s')
default_date_format = '%Y-%m-%d %H:%M:%S'

logging.basicConfig(format=default_logger_format, level=logging.INFO,
                    datefmt=default_date_format)


logger = logging.getLogger('genewalk')
# Remove info about missing package from gensim
logging.getLogger('gensim.summarization.textcleaner').setLevel(logging.WARNING)


__version__ = '1.2.2'
