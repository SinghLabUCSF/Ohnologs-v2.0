
# OHNOLOGS v2.0

##### OHNOLOGS v2.0 is a comprehensive repository for ohnologs retained from Whole Genome Duplications (WGDs) in the vertebrate genomes.    

> This repository contains all the codes and data used to construct [OHNOLOGS v2.0](http://ohnologs.curie.fr/) server.  

More information can be found in the directories  below:

* [0_Organisms_In_Ohnolog2.0](0_Organisms_In_Ohnolog2.0)
   Organisms and the outgroups included in OHNOLOGS v2.0.
   
* [1_All_Genes](1_All_Genes)
  Details on all the genes, Gene Ontology terms, chromosomal locations etc.
  
* [2_Paralogs](2_Paralogs)
  Details on candidate paralogs, their duplication timing, and consensus duplication timing from Ensembl.
  
* [3_Orthologs](3_Orthologs)
  Getting and preparing the orthologs from Ensembl among vertebrates and between vertebrates and outgroups for synteny analysis.
  
* [4_All_PC_Gene_Seqs](4_All_PC_Gene_Seqs)
  Including Amphioxus as an outgroup which is not covered by Ensembl.
  
* [5_Run_synteny](5_Run_synteny)
 Run the synteny analysis to identify all ohnologs and their confidence score (q-score).
 
* [6_CombineAllVertebrates](6_CombineAllVertebrates)
Combining ohnologs from all vertebrates by taking a phylogenetically biased averaging of q-scores.

* [7_FilterOhnologs](7_FilterOhnologs)
  Filter candidate ohnologs based on the three predefined (strict, intermediate and relaxed) q-score criteria.
  
* [8_GenerateOhnoFamilies](8_GenerateOhnoFamilies)
Generating ohnolog families for each criteria using ohnolog pairs.


> #### *If you use the data or code in the repository, please cite:*
> 
> * Singh P.P. and Isambert H. (2019) **OHNOLOGS v2: A comprehensive resource for the genes retained from whole genome duplication in vertebrates**. [bioRxiv](https://www.biorxiv.org/content/10.1101/717124v1).
> 
> * Singh P.P., Arora J. and Isambert H. (2015) **Identification of ohnolog genes originating from whole genome duplication in early vertebrates, based on synteny comparison across multiple genomes.** [PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004394).

