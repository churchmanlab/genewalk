from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def main():
    install_list = ['numpy', 'pandas', 'networkx>=2.1', 'gensim', 'goatools',
                    'indra>=1.14.1', 'scipy>=1.3.0']

    setup(name='genewalk',
          version='1.2.2',
          description='Determine gene function based on network embeddings.',
          long_description=long_description,
          long_description_content_type='text/markdown',
          author='Robert Ietswaart',
          author_email='robert_ietswaart@hms.harvard.edu',
          url='https://github.com/churchmanlab/genewalk',
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
          keywords=['gene function', 'network', 'embedding'],
          packages=find_packages(),
          install_requires=install_list,
          tests_require=['nose'],
          include_package_data=True,
          entry_points={'console_scripts': ['genewalk = genewalk.cli:main']},
        )


if __name__ == '__main__':
    main()
