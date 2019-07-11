# GeneWalk

GeneWalk determines for individual genes the functions that are relevant in a
particular experimental condition. GeneWalk quantifies the similarity between
vector representations of a gene and annotated GO terms through representation
learning with random walks on a condition-specific gene regulatory network.
Similarity significance is determined through comparison with node similarities
from randomized networks. 

## Install GeneWalk
To install the latest release of GeneWalk (preferred):
```
pip install genewalk
```
To install the latest code from Github (typically ahead of releases):
```
pip install git+https://github.com/churchmanlab/genewalk
```

## Using GeneWalk

### GeneWalk command line interface
Once installed, GeneWalk can be run from the command line as `genewalk`, with
a set of required and optional arguments. The required arguments include the
project name, a path to a text file containing a list of genes, and an argument
specifying the types of genes in the file.

Below is the full documentation of the command line interface:

```
genewalk [-h] --project PROJECT --genes GENES --id_type
              {hgnc_symbol,hgnc_id,mgi_id}
              [--stage {all,node_vectors,null_distribution,statistics}]
              [--base_folder BASE_FOLDER]
              [--network_source {pc,indra,edge_list,sif}]
              [--network_file NETWORK_FILE] [--nproc NPROC] [--nreps NREPS]
              [--alpha_fdr ALPHA_FDR] [--save_dw SAVE_DW]


required arguments:
  --project PROJECT     A name for the project which determines the folder
                        within the base folder in which the intermediate and
                        final results are written. Must contain only
                        characters that are valid in folder names.
  --genes GENES         Path to a text file with a list of differentially
                        expressed genes. Thetype of gene identifiers used in
                        the text file are provided in the id_type argument.
  --id_type {hgnc_symbol,hgnc_id,mgi_id}
                        The type of gene IDs provided in the text file in the
                        genes argument. Possible values are: hgnc_symbol,
                        hgnc_id, and mgi_id.

optional arguments:
  --stage {all,node_vectors,null_distribution,statistics}
                        The stage of processing to run.
  --base_folder BASE_FOLDER
                        The base folder used to store GeneWalk temporary and
                        result files for a given project.
  --network_source {pc,indra,edge_list,sif}
                        The source of the network to be used.Possible values
                        are: pc, indra, edge_list, and sif. In case of indra,
                        edge_list, and sif, the network_file argument must be
                        specified.
  --network_file NETWORK_FILE
                        If network_source is indra, this argument points to a
                        Python pickle file in which a list of INDRA Statements
                        constituting the network is contained. In case
                        network_source is edge_list or sif, the network_file
                        argument points to a text file representing the
                        network.
  --nproc NPROC         The number of processors to use in a multiprocessing
                        environment.
  --nreps_graph NREPS_GRAPH
                        The number of repeats to run when calculating node
                        vectors on the "real" network graph.
  --nreps_null NREPS_NULL
                        The number of repeats to run when calculating node
                        vectors on the random network graphs for constructing
                        the null distribution.
  --alpha_fdr ALPHA_FDR
                        The false discovery rate to use when calculating the
                        final statistics.
  --save_dw SAVE_DW     If True, the full DeepWalk object for each repeat is
                        saved in the project folder. This can be useful for
                        debugging but the files are typically very large.
  --random_seed RANDOM_SEED
                        If given, the random number generator will be seeded
                        with the given value. This should only be used if the
                        goal is to deterministically reproduce a prior result
                        obtained with the same random seed.

```

Example
```bash
genewalk --project test1 --genes gene_list.txt --id_type hgnc_symbol
```

### Gene list file
GeneWalk always requires a text file with a list of genes as input. These are
typically differentially expressed genes in a given experimental setting.
GeneWalk supports gene list files containing HGNC symbols, HGNC IDs, or MGI
IDs. Each line in the file contains a gene identifier of one of these types.


### Output files
GeneWalk automatically creates a `genewalk` folder in the user's home folder.
When running GeneWalk, one of the required inputs is a project name.
A sub-folder is created for the given project name where all intermediate and
final results are stored. The files stored in the project folder are:
- genewalk_results.csv - The main results table, a comma-separated values text file.
- genes.pkl - A processed representation of the given gene list.
- multi_graph.pkl - A networkx MultiGraph which was assembled based on the
given list of genes, an interaction network, GO annotations, and the GO
ontology.
- deep_walk_*.pkl - A DeepWalk object for each analysis repeat on the graph.
- deep_walk_rand_*.pkl - A DeepWalk object for each analysis repeat on a random graph.
- deep_walk_node_vectors_*.pkl - A set of learned node vectors for each analysis repeat for the graph.
- deep_walk_node_vectors_rand_*.pkl - A set of learned node vectors for each analysis repeat for a random graph.
- deep_walk_rand_simdists.pkl - Distributions constructed from repeats.

### Stages of analysis
Given a list of genes, GeneWalk runs three stages of analysis:
1. Assembling a GeneWalk network and learning node vector representations
by running DeepWalk on this network, for a specified number of repeats.
2. Learning random node vector representations by running DeepWalk on a set of
randomly scrambled versions of the GeneWalk network, for a specified number of
repeats.
3. Calculating statistics of similarities between genes and GO terms, and
outputting  the GeneWalk results in a table.

GeneWalk can either be run once to complete all these stages, or called separately
for each stage.
