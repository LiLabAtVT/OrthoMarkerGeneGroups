# Convert species marker genes into OMGs by OrthoFinder

### Install:
Install OrthoFinder on Mac by using biconda: <br>
`conda install orthofinder`

### Run OrthoFinder:
Prepare input files in FASTA format for each species in the fasta folder
Run OrthoFinder on the input fasta files: <br>
`OrthoFinder/orthofinder -f OrthoFinder/fasta`

### Result:
The output data is in the zip Orthogroups.tsv.zip file.

### Note:
Make sure the gene IDs in your marker gene list match the gene IDs used in the FASTA headers given to OrthoFinder. Mismatched or differently formatted IDs (e.g. with version suffixes, transcript vs. gene IDs, or different naming conventions) will prevent your marker genes from being mapped to their OMGs.

For example, if your marker gene list uses `AT1G01010` but the FASTA header is written as `AT1G01010.1` (transcript suffix) or `>AT1G01010.1 | Symbol: NAC001`, the IDs will not match. Standardize both so they use the same form. Examples for a few species:

```
Species                  Marker gene list      FASTA header
Arabidopsis thaliana     AT1G01010             >AT1G01010
Solanum lycopersicum     Solyc12g009745        >Solyc12g009745
Populus trichocarpa      Potri_T125208         >Potri_T125208
```

A common mismatch is a trailing transcript suffix, e.g. marker `Solyc12g009745` vs. header `>Solyc12g009745.2.1` , strip the suffix so both sides are identical.

### Reference:
Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20, 1-14.
