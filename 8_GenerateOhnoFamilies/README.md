## Generate Ohnolog families using depth first search

These scripts generate ohnolog families by combining pairs for all the organisms. 

* `0_create_dirs.pl`: Create directpries for all organisms.

* `1_DepthFirstSearchOhnolgFamilies_cl.pl`: Generate ohnolog families using depth first search. It will use all the pairs and get all ohnologs connected to these pairs until no new ohnolog memebers are found. 

* `1b_filterSize1families_cl.pl`: Filter our families if they are all merged together. This can happen if they lie within the smallest window after combining all vertebrates.

* `2_filterAncientDuplicates_cl.pl`: Filter old small scale duplicates (SSDs) that occurred before WGD. Ancient SSDs are separated by a pipe.

* `3_filterRecentDuplicates_cl.pl`: Filter recent SSDs that occurred before WGD. Recent SSDs are separated by comma.

* `3b_getCliqueFamilies_cl.pl`: Getting families that for clique i.e. each memeber is connected by every other member.

* `4_getGeneSymols.pl`: Convert Ensembl Ids to symbols.

* `run_all.pl`: Run the piepeline for all or selected organisms.

> Note that we remove SSDs age information from the OHNOLOGS sever for simplicity. On the server all SSDs within an ohnon family are separated by a comma regardless of age.

