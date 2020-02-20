# GeneWalk

[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Documentation](https://readthedocs.org/projects/genewalk/badge/?version=latest)](https://genewalk.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/genewalk.svg)](https://badge.fury.io/py/genewalk)
[![Python 3.5+](https://img.shields.io/pypi/pyversions/genewalk.svg)](https://www.python.org/downloads/release/python-357/)

GeneWalk determines for individual genes the functions that are relevant in a
particular biological context and experimental condition. GeneWalk quantifies
the similarity between vector representations of a gene and annotated GO terms
through representation learning with random walks on a condition-specific gene
regulatory network. Similarity significance is determined through comparison
with node similarities from randomized networks.

## Install GeneWalk
To install the latest release of GeneWalk (preferred):
```
pip install genewalk
```
To install the latest code from Github (typically ahead of releases):
```
pip install git+https://github.com/churchmanlab/genewalk.git
```

## Using GeneWalk

### Gene list file
GeneWalk always requires as input a text file containing a list with genes of interest
relevant to the biological context. For example, differentially expressed genes
from a sequencing experiment that compares an experimental versus control condition.
GeneWalk supports gene list files containing HGNC human gene symbols,
HGNC IDs, human Ensembl gene IDs, or MGI mouse gene IDs. Each line in the file
contains a gene identifier of one of these types.

### GeneWalk command line interface
Once installed, GeneWalk can be run from the command line as `genewalk`, with
a set of required and optional arguments. The required arguments include the
project name, a path to a text file containing a list of genes, and an argument
specifying the types of genes in the file.

Example
```bash
genewalk --project context1 --genes gene_list.txt --id_type hgnc_symbol
```

Below is the full documentation of the command line interface:

```
genewalk [-h] [--version] --project PROJECT --genes GENES --id_type
              {hgnc_symbol,hgnc_id,mgi_id,ensembl_id}
              [--stage {all,node_vectors,null_distribution,statistics}]
              [--base_folder BASE_FOLDER]
              [--network_source {pc,indra,edge_list,sif}]
              [--network_file NETWORK_FILE] [--nproc NPROC] [--nreps NREPS]
              [--alpha_fdr ALPHA_FDR] [--save_dw SAVE_DW]
              [--random_seed RANDOM_SEED]


required arguments:
  --version             Print the version of GeneWalk and exit.
  --project PROJECT     A name for the project which determines the folder
                        within the base folder in which the intermediate and
                        final results are written. Must contain only
                        characters that are valid in folder names.
  --genes GENES         Path to a text file with a list of differentially
                        expressed genes. Thetype of gene identifiers used in
                        the text file are provided in the id_type argument.
  --id_type {hgnc_symbol,hgnc_id,ensembl_id,mgi_id,entrez_human,entrez_mouse}
                        The type of gene IDs provided in the text file in the
                        genes argument. Possible values are: hgnc_symbol,
                        hgnc_id, ensembl_id, mgi_id, entrez_human and
                        entrez_mouse.

optional arguments:
  --stage {all,node_vectors,null_distribution,statistics}
                        The stage of processing to run. Default: all
  --base_folder BASE_FOLDER
                        The base folder used to store GeneWalk temporary and
                        result files for a given project. Default:
                        ~/genewalk
  --network_source {pc,indra,edge_list,sif}
                        The source of the network to be used.Possible values
                        are: pc, indra, edge_list, and sif. In case of indra,
                        edge_list, and sif, the network_file argument must be
                        specified. Default: pc
  --network_file NETWORK_FILE
                        If network_source is indra, this argument points to a
                        Python pickle file in which a list of INDRA Statements
                        constituting the network is contained. In case
                        network_source is edge_list or sif, the network_file
                        argument points to a text file representing the
                        network.
  --nproc NPROC         The number of processors to use in a multiprocessing
                        environment. Default: 1
  --nreps_graph NREPS_GRAPH
                        The number of repeats to run when calculating node
                        vectors on the GeneWalk graph. Default: 3
  --nreps_null NREPS_NULL
                        The number of repeats to run when calculating node
                        vectors on the random network graphs for constructing
                        the null distribution. Default: 3
  --alpha_fdr ALPHA_FDR
                        The false discovery rate to use when outputting the
                        final statistics table. If 1 (default), all
                        similarities are output, otherwise only the ones whose
                        false discovery rate are below this parameter are
                        included. Default: 1
  --save_dw SAVE_DW     If True, the full DeepWalk object for each repeat is
                        saved in the project folder. This can be useful for
                        debugging but the files are typically very large.
                        Default: False
  --random_seed RANDOM_SEED
                        If provided, the random number generator is seeded
                        with the given value. This should only be used if the
                        goal is to deterministically reproduce a prior result
                        obtained with the same random seed.

```


### Output files
GeneWalk automatically creates a `genewalk` folder in the user's home folder
(or the user specified base_folder).
When running GeneWalk, one of the required inputs is a project name.
A sub-folder is created for the given project name where all intermediate and
final results are stored. The files stored in the project folder are:
- **`genewalk_results.csv`** - The main results table, a comma-separated values text file. See below for detailed description.
- `genes.pkl` - A processed representation of the given gene list, in Python pickle (.pkl) binary file format.
- `multi_graph.pkl` - A networkx MultiGraph resembling the GeneWalk network which was assembled based on the
given list of genes, an interaction network, GO annotations, and the GO ontology.
- `deepwalk_node_vectors_*.pkl` - A set of learned node vectors for each analysis repeat for the graph.
- `deepwalk_node_vectors_rand_*.pkl` - A set of learned node vectors for each analysis repeat for a random graph.
- `genewalk_rand_simdists.pkl` - Distributions constructed from repeats.
- `deepwalk_*.pkl` - A DeepWalk object for each analysis repeat on the graph
(only present if save_dw argument is set to True).
- `deepwalk_rand_*.pkl` - A DeepWalk object for each analysis repeat on a random graph
(only present if save_dw argument is set to True).


### GeneWalk results file description
`genewalk_results.csv` is the main GeneWalk output table, a comma-separated values text file
with the following column headers:
- hgnc_id - human gene HGNC identifier.
- **hgnc_symbol** - human gene symbol.
- **go_name** - GO term name.
- go_id - GO term identifier.
- go_domain - Ontology domain that GO term belongs to
(biological process, cellular component or molecular function).
- ncon_gene - number of connection to gene in GeneWalk network.
- ncon_go - number of connections to GO term in GeneWalk network.
- **mean_padj** - mean false discovery rate (FDR) adjusted p-value of the similarity between gene and GO term.
This is the key statistic indicating how relevant the GO term (function) is to the gene in the
particular biological context or tested condition. GeneWalk determines an adjusted p-value with
Benjamini Hochberg FDR correction for multiple tested of all connected GO term for each
nreps_graph repeat analysis. The value presented here is the average over all p-adjust values
from each repeat analysis. 
- cilow_padj - lower bound of 95% confidence interval on mean_padj estimate from the nreps_graph repeat analyses.
- ciupp_padj - upper bound of 95% confidence interval on mean_padj estimate.
- mean_pval - mean p-values of gene - GO term similarities, not FDR corrected for multiple testing.
- cilow_pval - lower bound of 95% confidence interval on mean_pval estimate.
- ciupp_pval - upper bound of 95% confidence interval on mean_pval estimate.
- mean_sim - mean of gene - GO term similarities.
- sem_sim - standard error on mean_sim estimate.
- mgi_id, ensembl_id, mgi_id, entrez_human or entrez_mouse - in case one of
  these gene identifiers were provided as input, the GeneWalk results table
  starts with an additional column to indicate the gene identifiers. In the
  case of mouse genes, the corresponding hgnc_id and hgnc_symbol resemble its
  human ortholog gene used for the GeneWalk analysis.


### Run time and stages of GeneWalk algorithm
Recommended number of processors (optional argument: nproc) for a short (1-2h)
run time is 4:
```bash
genewalk --project context1 --genes gene_list.txt --id_type hgnc_symbol --nproc 4
```
By default GeneWalk will run with 1 processor, resulting in a longer overall
run time: 6-12h.
Given a list of genes, GeneWalk runs three stages of analysis:
1. Assembling a GeneWalk network and learning node vector representations
by running DeepWalk on this network, for a specified number of repeats.
Typical run time: one to a few hours.
2. Learning random node vector representations by running DeepWalk on a set of
randomized versions of the GeneWalk network, for a specified number of
repeats. Typical run time: one to a few hours.
3. Calculating statistics of similarities between genes and GO terms, and
outputting  the GeneWalk results in a table. Typical run time: a few minutes.

GeneWalk can either be run once to complete all these stages (default), or
called separately for each stage (optional argument: stage).  Recommended
memory availability on your operating system: 16Gb or 32Gb RAM.  GeneWalk
outputs the uncertainty (95% confidence intervals) of the similarity
significance (mean p-adjust). Depending on the context-specific network
topology, this uncertainty can be large for individual gene - function
associations. However, if overall the uncertainties turn out very large, one
can set the optional arguments nreps_graph to 10 (or more) and nreps_null to 10
to increase the algorithm's precision. This comes at the cost of an increased
run time.


### Further documentation
For a tutorial and more general information see the
[GeneWalk website](http://churchman.med.harvard.edu/genewalk).
For further code documentation see our [readthedocs page](https://genewalk.readthedocs.io).


### Citation
Robert Ietswaart, Benjamin M. Gyori, John A. Bachman, Peter K. Sorger, and
L. Stirling Churchman
*GeneWalk identifies relevant gene functions for a biological context using network
representation learning* (2019), [BioRxiv; 755579](https://www.biorxiv.org/content/10.1101/755579v2).


### Funding
This work was supported by National Institutes of Health grant 5R01HG007173-07
(L.S.C.), EMBO fellowship ALTF 2016-422 (R.I.), and DARPA grants W911NF-15-1-0544
and W911NF018-1-0124 (P.K.S.).
