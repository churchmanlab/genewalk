# GeneWalk

`source genewalkenv/bin/activate`
`outDir=~/genewalk`
`mkdir -p ${outDir}/LogErr`
`python get_indra_stmts.py` #to BG/JB: this line needs more specification
`python -u get_node_vectors.py --path ${outDir} > ${outDir}/LogErr/get_node_vectors.log`
`python -u get_null_distributions.py --path ${outDir} > ${outDir}/LogErr/get_null_distributions.log`
