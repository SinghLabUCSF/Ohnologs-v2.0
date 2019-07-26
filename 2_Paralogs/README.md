## Download and filter paralogs

1. Use  `1_get_biomaRt_Paralogs.R`  or `1.1_get_biomaRt_Paralogs_FileTypes.R` to get paralogs for all genes types or for each gene type one by one for large datastes.

2.  `2_getAllNodes.pl`  and `3_getReconsiledNodes.pl` combine multiple versions from v80-86, to get the consensus duplication time which is more accurate.

3. `3b_filter_multi-copy_paralogs.pl` filters the non-coding RNA genes that have more than 30 paralogs. This is to limit the number of paralogs  and improves computation time, which otherwise is impractical. 

4. Filter paralogs duplicated at the base of vertebrates or fish for 2R and 2R/3R WGS respectively using   `4_filterReconsiledNodes.pl` and `4b_filterReconsiledNodes_Fish.pl`.

> Candidate duplication nodes used:
> * For 2R WGD: Vertebrata, Euteleostomi, Chordata, Sarcopterygii or Neopterygii.    
> * For 3R-WGD: FishWGD, Clupeocephala, Acanthomorphata. FishWGD has paralog candidates for the consensus node that is duplicated at the time of 3R-WGD.
