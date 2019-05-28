# GeneWalk

### Install GeneWalk
To install the latest release of GeneWalk, do
```
pip install genewalk
```
To install the latest code from Github (typically ahead of releases):
```
pip install git+https://github.com/churchmanlab/genewalk
```
Optionally, to use INDRA's sources to collect information
and assemble a network, do
```
pip install indra
```

Note that it is recommended to use Python virtual environments when installing
GeneWalk. You can start a new environment before installation as
```
virtualenv genewalkenv --system-site-packages
source genewalkenv/bin/activate
```

### Set output directory
`outDir=~/genewalk`  
`mkdir -p ${outDir}/LogErr`  

### Run GeneWalk
`python get_indra_stmts.py` #################to BG/JB: this line needs more specification / alteration  
`python -u get_node_vectors.py --path ${outDir} > ${outDir}/LogErr/get_node_vectors.log`  
`python -u get_null_distributions.py --path ${outDir} > ${outDir}/LogErr/get_null_distributions.log`  
`python -u genewalk.py --path ${outDir} > ${outDir}/LogErr/genewalk.log`
