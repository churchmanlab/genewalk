"""This script creates a SIF file for a specific set of genes, based on a
pre-assembled network of interactions produced by INDRA. It loads the INDRA
interactions as a pandas DataFrame and, filters it to the genes of interest,
and adds some additional information before exporting as SIF."""
import pandas
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Choose a file with a list of genes to get a SIF for.')
    parser.add_argument('--genes')
    args = parser.parse_args()
    genes = []
    with open(args.genes, 'r') as fh:
        # Get the HGNC IDs from the list of genes file, assuming that
        # each line looks like HGNC:123
        genes = [l.strip().split(':')[1] for l in fh.readlines()]
