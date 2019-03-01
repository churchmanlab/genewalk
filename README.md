# GeneWalk

### Initialize new genewalk virtual environment
`virtualenv genewalkenv --system-site-packages`  
`source genewalkenv/bin/activate`

### Install dependencies
`pip install pandas`  
`pip install --upgrade numpy`  
`pip install --upgrade scipy`  
`pip install git+https://github.com/sorgerlab/indra.git`  
`pip install goatools`  
`pip install --upgrade gensim`  
`pip install statsmodels`  
`pip install genewalk` #######################UNFINISHED, to BG/JB: please help adjust so that this could work 

### Set output directory
`outDir=~/genewalk`  
`mkdir -p ${outDir}/LogErr`  
`mkdir -p ${outDir}/GO`

### Download GO ontology and annotations
`cd ${outDir}/GO`

Linux:  `wget http://snapshot.geneontology.org/ontology/go.obo`  
`wget http://geneontology.org/gene-associations/goa_human.gaf.gz`

or with Mac OS:  `curl -O http://snapshot.geneontology.org/ontology/go.obo`  
`curl -O http://geneontology.org/gene-associations/goa_human.gaf.gz`

`gunzip http://geneontology.org/gene-associations/goa_human.gaf.gz`

### Run GeneWalk
`python get_indra_stmts.py` #################to BG/JB: this line needs more specification / alteration  
`python -u get_node_vectors.py --path ${outDir} > ${outDir}/LogErr/get_node_vectors.log`  
`python -u get_null_distributions.py --path ${outDir} > ${outDir}/LogErr/get_null_distributions.log`  
`python -u genewalk.py --path ${outDir} > ${outDir}/LogErr/genewalk.log`