# Convert species marker genes into OMGs by OrthoFinder

### Install:
Install OrthoFinder on Mac by using biconda: <br>
`conda install orthofinder`

### Run OrthoFinder:
Prepare input files in FASTA format for each species in the fasta folder
Run OrthoFinder on the input fasta files: <br>
`OrthoFinder/orthofinder -f OrthoFinder/fasta`

### Result:
The output data is in the zip Orthogroups_091023_cleaned.tsv.zip file.

### Reference:
Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20, 1-14.
