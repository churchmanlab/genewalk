# GeneWalk

Optional: it is recommended to use Python virtual environments when installing
GeneWalk. You can start a new environment before installation as
```
virtualenv genewalkenv --system-site-packages
source genewalkenv/bin/activate
```

### Install GeneWalk
To install the latest release of GeneWalk (preferred), do
```
pip install genewalk
```
To install the latest code from Github (typically ahead of releases):
```
pip install git+https://github.com/churchmanlab/genewalk
```



### Set output directory
`outDir=~/genewalk`  
`mkdir -p ${outDir}/LogErr`  

### Run GeneWalk
with the Pathway Commons data source (currently default):
`python -u get_node_vectors.py --path ${outDir} > ${outDir}/LogErr/get_node_vectors.log`  
`python -u get_null_distributions.py --path ${outDir} > ${outDir}/LogErr/get_null_distributions.log`  
`python -u genewalk.py --path ${outDir} > ${outDir}/LogErr/genewalk.log`

### Optional: GeneWalk with INDRA data source
to use INDRA's sources to collect information and assemble a network, do
```
pip install indra
```

`python get_indra_stmts.py`
`python -u get_node_vectors.py --path ${outDir} > ${outDir}/LogErr/get_node_vectors.log`  
`python -u get_null_distributions.py --path ${outDir} > ${outDir}/LogErr/get_null_distributions.log`  
`python -u genewalk.py --path ${outDir} > ${outDir}/LogErr/genewalk.log`
