import logging

logging.basicConfig(format=('%(levelname)s: [%(asctime)s] %(name)s'
                            ' - %(message)s'),
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')#,
#                     handlers=logging.FileHandler('/n/groups/churchman/ri23/test.log'))#somehow this does not work yet

__version__ = '0.0.1'
