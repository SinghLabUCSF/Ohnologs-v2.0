
## Combine ohnologs from all windows for self and outgroup comparison

* Window sizes used are 100, 200, 300, 400 and 500; for both self and outgroup comparisons.
* Total 29 organisms for 2R-WGD, and 4 fish for 3R-WGD (>1000 total comparisons)
* `0_create_dirs.pl`: Create a directory for each organism to hold results. Organism name are in the script.
* `1_CombinePairsFromDifferentWindows_Self_All.pl`: Combine all windows for self and filter 2-way relationships.
* `2_CombinePairsFromDifferentOutgroupWindow_All.pl`: Combine all windows for each outgroup.
* `3_CombinePairsFromDifferentOutgroup_All.pl`: Combine all outgroup pairs from the step above.