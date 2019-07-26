## Download orthologs 

1. `1_get_biomaRt_Orthologs.R` downloads orthologs between outgroups and WGD organisms for 2R and 3R WGDs using the organism lists. 

2. Use `1b_get_biomaRt_Orthologs-among-ALL-vertebrates.R`
to get all orthologs. 

3. Add Amphioxus using `2_getAmphioxusOrthologs.pl` from the BLASTp results dirs.
   
4. Filter the Ensmebl orthologs to remove blank lines and NAs using `3_filterEnsemblOrthologs.pl`.   
