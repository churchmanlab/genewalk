# GeneWalk

[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Documentation](https://readthedocs.org/projects/genewalk/badge/?version=latest)](https://genewalk.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/genewalk.svg)](https://badge.fury.io/py/genewalk)
[![Python 3.6+](https://img.shields.io/pypi/pyversions/genewalk.svg)](https://www.python.org/downloads/release)

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

GeneWalk uses a number of resource files that it downloads as needed during
runtime. To optionally pre-download these resource files in the default resource folder,
the command
```
python -m genewalk.resources
```
can be run.

## Using GeneWalk

### Gene list file
GeneWalk always requires as input a text file containing a list with genes of interest
relevant to the biological context. For example, differentially expressed genes
from a sequencing experiment that compares an experimental versus control condition.
GeneWalk supports gene list files containing HGNC human gene symbols,
HGNC IDs, human Ensembl gene IDs, MGI mouse gene IDs, RGD rat gene IDs, or
human or mouse entrez IDs. Each line in the file contains a gene identifier
of one of these types.

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
              {hgnc_symbol,hgnc_id,ensembl_id,mgi_id,rgd_id,entrez_human,entrez_mouse}
              [--stage {all,node_vectors,null_distribution,statistics}]
              [--base_folder BASE_FOLDER]
              [--network_source {pc,indra,edge_list,sif,sif_annot,sif_full}]
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
  --id_type {hgnc_symbol,hgnc_id,ensembl_id,mgi_id,rgd_id,entrez_human,entrez_mouse}
                        The type of gene IDs provided in the text file in the
                        genes argument. Possible values are: hgnc_symbol,
                        hgnc_id, ensembl_id, mgi_id, rgd_id, entrez_human and
                        entrez_mouse.

optional arguments:
  --stage {all,node_vectors,null_distribution,statistics,visual}
                        The stage of processing to run. Default: all
  --base_folder BASE_FOLDER
                        The base folder used to store GeneWalk temporary and
                        result files for a given project. Default:
                        ~/genewalk
  --network_source {pc,indra,edge_list,sif,sif_annot,sif_full}
                        The source of the network to be used.Possible values
                        are: pc, indra, edge_list, sif, sif_annot, and
                        sif_full. In case of indra, edge_list, sif, sif_annot,
                        and sif_full, the network_file argument must be
                        specified. Default: pc
  --network_file NETWORK_FILE
                        If network_source is indra, this argument points to a
                        Python pickle file in which a list of INDRA Statements
                        constituting the network is contained. In case
                        network_source is edge_list, sif, sif_annot, or
                        sif_full, the network_file argument points to a text
                        file representing the network.
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
                        For visualization a default value of 0.1 for both global
                        and gene-specific plots is used. Lower this value to 
                        increase the stringency of the regulator gene selection 
                        procedure.
  --dim_rep DIM_REP     Dimension of vector representations (embeddings). This 
                        value should only be increased if genewalk with the 
                        default value generates no statistically significant 
                        results, for instance with very large (>2500) input 
                        gene lists. Alternatively, it can be decreased in case 
                        (nearly) all GO annotations are significant, for 
                        instance with very short gene lists. Default: 8
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

### Figure files
GeneWalk also automatically generates figures to visualize its results in the
project/figures sub-folder:
- **`index.html`**: an HTML page that includes all the figures generated, as
  described below.
- barplots with GO annotations ranked by relevance for each input gene that
  GeneWalk was able to generate results for. The filenames contain the
  corresponding human gene symbol and input gene id: `barplot_[symbol]_[gene
  id]_x_mlog10global_padj_y_GO.png`.
- `regulators_x_gene_con_y_frac_rel_go(.png and .pdf)`: scatter plot to
  identify regulator genes of interest. These have a large gene connectivity
  and high fraction of relevant GO annotations. For more information see our
  publication.
- `genewalk_regulators.csv`: list with regulator genes that are named in the
  regulators scatterplot.
- `moonlighters_x_go_con_y_frac_rel_go(.png and .pdf)`: scatter plot to
  identify moonlighting genes: genes with many GO annotations of which a low
  fraction are relevant. For more information see our publication.
- `genewalk_moonlighters.csv`: list with moonlighting genes that are named in
  the moonlighting scatterplot.
- `genewalk_scatterplots.csv`: data corresponding to the regulator and
  moonlighter scatter plots.  This file can be used for further gene
  prioritization analyses.


### GeneWalk results file description
`genewalk_results.csv` is the main GeneWalk output table, a comma-separated values text file
with the following column headers:
- hgnc_id - human gene HGNC identifier.
- **hgnc_symbol** - human gene symbol.
- **go_name** - GO term name.
- go_id - GO term identifier.
- go_domain - Ontology domain that GO term belongs to
(biological process, cellular component or molecular function).
- ncon_gene - number of connections to gene in GeneWalk network.
- ncon_go - number of connections to GO term in GeneWalk network.
- **global_padj** - false discovery rate (FDR) adjusted p-value of the 
similarity between gene and GO term, when correcting for testing over all 
gene-GO term pairs present in the output file.
This is the key statistic that indicates how relevant the gene-GO term pair 
(gene function) is in the particular biological context or tested condition. 
Global_padj should be used for global analyses that
consider all the GeneWalk output simultaneously, such as gene prioritization
procedures. GeneWalk determines an adjusted p-value with Benjamini Hochberg FDR 
correction for multiple testing of all connected GO term for each 
nreps_graph repeat analysis. The value presented here is the average (mean 
estimate) over all p-adjust values from all nreps_graph repeat analyses. 
- **gene_padj** - FDR adjusted p-value of the similarity between gene and 
GO term, when correcting for multiple testing over all GO annotations of 
that gene. This the key statistic when investigating the functions of one 
(or a few) pre-defined gene(s) of interest. Gene_padj determines the statistical 
significance of each GO annotation (function) and gene_padj can be used to 
sensitively rank GO annotations to reflect the relevance to the gene of interest
in the particular biological context or tested condition. When you consider all
(or many) input genes simultaneously, use global_padj instead. Average 
over nreps_graph repeat runs as for global_padj. 
- pval - p-value of gene - GO term similarity, not corrected for multiple
hypothesis testing. Average over nreps_graph repeat runs.
- sim - gene - GO term (cosine) similarity, average over nreps_graph repeat runs.
- sem_sim - standard error on sim (mean estimate).
- cilow_global_padj - lower bound of 95% confidence interval on global_padj 
(mean estimate) from the nreps_graph repeat analyses.
- ciupp_global_padj - upper bound of 95% confidence interval on global_padj.
- cilow_gene_padj - lower bound of 95% confidence interval on gene_padj
(mean estimate) from the nreps_graph repeat analyses.
- ciupp_gene_padj - upper bound of 95% confidence interval on gene_padj.
- cilow_pval - lower bound of 95% confidence interval on pval (mean estimate)
from the nreps_graph repeat analyses.
- ciupp_pval - upper bound of 95% confidence interval on pval.
- mgi_id, rgd_id, ensembl_id, entrez_human or entrez_mouse - in case one of
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
4. Visualization of the GeneWalk results generated in the project/figures subfolder.
Typical run time: 1-10 mins depending on the number of input genes.

GeneWalk can either be run once to complete all these stages (default), or
called separately for each stage (optional argument: stage).  Recommended
memory availability on your operating system: 16Gb or 32Gb RAM.  GeneWalk
outputs the uncertainty (95% confidence intervals) of the similarity
significance (global and gene p-adjust). Depending on the context-specific network
topology, this uncertainty can be large for individual gene - function
associations. However, if overall the uncertainties turn out very large, one
can set the optional arguments nreps_graph to 10 (or more) and nreps_null to 10
to increase the algorithm's precision. This comes at the cost of an increased
run time.


### Custom input networks
By default, GeneWalk uses the PathwayCommons network (`--network_source pc`)
to create a human gene network. It then automatically adds edges
reprenting GO annotations for input genes as well as relations between
GO terms. However, there are also options for supplying a custom network as
input, as follows.

The `--network_source sif/sif_annot/sif_full` options require supplying the
path to a simple interaction file (SIF) as the `--network_file` argument. Each
row of the SIF file consists of three comma-separated entries representing
source, relation type, and target. Genes in the SIF are assumed to be
human gene symbols (e.g., KRAS), and GO terms in the SIF (only in the
`sif_annot` and `sif_full` modes) use the GO: prefix (e.g., GO:0000186).
The relation type is not explicitly used by GeneWalk, and can be set
to an arbitrary label.

The difference between the `sif`,
`sif_annot`, and `sif_full` options are as follows:
- `sif`: In this case, the input SIF is assumed to contain only gene-gene
   relations. GO annotations for genes, as well as relations between
   GO terms are added automatically by GeneWalk.
- `sif_annot`: In this case, the input SIF is assumed to contain both
  gene-gene relations, and GO annotations for genes (i.e., rows where the
  source is a gene, and the target is a GO term). Relations between
  GO terms are then added automatically by GeneWalk.
- `sif_full`: In this case, the input SIF is assumed to contain all
  relations including gene-gene relations, GO annotations for genes,
  and relations between GO terms. GeneWalk doesn't add any further
  edges.

The `--network_source edge_list` option is a simplified version of the
`--network_source sif` option. It requires a `--network_file` argument
which contains rows with two columns each, a source and a target
(i.e., it omits the relation type column from SIF).

The `--network_source indra` option requires supplying the path to a Python
pickle file containing a list of INDRA Statements as the `--network_file`
argument. These statements can represent gene-gene, as well as gene-GO
relations from which network edges are derived. GO annotations, as well as
relations between GO terms are then added automatically by GeneWalk during
network construction.


### Further documentation
For a tutorial and more general information see the
[GeneWalk website](http://churchman.med.harvard.edu/genewalk).  
For further code documentation see our [readthedocs page](https://genewalk.readthedocs.io).


### Citation
Robert Ietswaart, Benjamin M. Gyori, John A. Bachman, Peter K. Sorger, and
L. Stirling Churchman  
*GeneWalk identifies relevant gene functions for a biological context using network
representation learning*,  
Genome Biology **22**, 55 (2021). [https://doi.org/10.1186/s13059-021-02264-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02264-8)  


### Funding
This work was supported by National Institutes of Health grant 5R01HG007173-07
(L.S.C.), EMBO fellowship ALTF 2016-422 (R.I.), and DARPA grants W911NF-15-1-0544
and W911NF018-1-0124 (P.K.S.).
