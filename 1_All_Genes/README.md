## Download genes and attribute lists from Ensembl 

Use biomaRt to download protein coding genes for each of the genome to be analyzed.

* Script: `1_get_biomaRt_GeneAttribute_files`.

* [species.txt](./species.txt) has the name of species list.

* Ensembl version used is v84

* Similarly Gene Ontology attributes have been downloaded using `2_get_biomaRt_GO_files.R` (see below for more details).

#### Prepare a sorted gene list file for the OHNOLOGS runs.

* In the earlier version of OHNOLOGS I removed genes not on assembled chromosomes ( from Ensembl karyotype pages for the 6 organisms). In v2, I keep all the chromosomes and scaffolds as they can have important genes that will otherwise be ignored.
  
* I remove mitochondrial gene sequences from all the organisms.
  
* I curated a Karyotype page in Ensembl for all the organisms to decide which genes to discard from the analysis. This information can be found in [ChromosomeCuration_Ens84.xlsx](ChromosomeCuration_Ens84.xlsx).
  
* Not all gene types have ortholog and paralog information. So I only keep following gene types:
  * miRNA or mirna
  * misc_RNA
  * protein_coding
  * rRNA
  * snRNA
  * snoRNA
  
* Use `3_processAllPCfile_with-scaffolds.pl` to prepare lists of all genes for OHNOLOGS.

* Amphioxus is processed separately. See [README](../4_All_PC_Gene_Seqs/README.md).
    
