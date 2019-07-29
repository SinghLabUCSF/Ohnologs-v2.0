## Combine all vertebrates 

> Here, we combine ohnologs from all vertebrates by taking a **phylogenetically biased q-score** average. 

1. Create directories for each organism `0_create_dirs.pl` 

2. `1_Combine_OhnoPairs_AllVertebrates_fish_3R.pl` and `1_Combine_OhnoPairs_AllVertebrates_tetra-fish_2R.pl`  to combine all vertebrate q-scores for 3R or 2R-WGD.
   
3. `2_getGeoMean_fish_3R_weighted-mean.pl` and `2_getGeoMean_tetra-fish_2R_weighted-mean.pl` calculate the phylogenetically weighted mean of q-scores across all vertebrates based on weights in (weights)[weights.txt].
   
4. `3_filterWithinWindows.pl` filters the genes that may be combined due to their ortholog in other vertebrates but that lie within the smallest window of 100 genes.

5. `4_addMissingInfo.pl` adds the duplication time and gene names if it's missing for some genes. 	