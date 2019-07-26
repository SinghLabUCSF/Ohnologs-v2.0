# Set wd to current file location
setwd("/Volumes/Mybook_2/Ohnologs/Synteny_All_2016_03_09/2_Paralogs")
# read organism names
org = readLines("species.txt")
#org = readLines("speciestest.txt")
versions = readLines("ensembl_versions_80-86.txt")

# load biomart libraray
library("biomaRt")


# This is to test and print some of the importnat variables
# This is for test purpose only, and we can comment it once we are sure about the data
# ------------------------------------------------------------------------------------
#listMarts()                  # List the available BioMarts and the version
#ensembl=useMart("ensembl")   # select to query the Ensembl BioMart database
#listDatasets(ensembl)        # for Ensembl every species is a different dataset. List all available datasets
# Write all datasets in a file
#write.table(listDatasets(ensembl), file = "ensembldatasets.txt", quote = F, sep = "\t", row.names = F)
#humanensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl) # select a datset
#listAttributes(humanensembl) # List attributes for the selected dataset
#write.table(listAttributes(humanensembl), file = "attributes.txt", quote = F, sep = "\t", row.names = F)
#  rm(list=ls()) # One we have seen everything, clear the workspace before I run the dataset
# ------------------------------------------------------------------------------------
#
# Ensembl data can be downloaded using getBM function, getBM has three arguments that 
# need to be introduced: filters, attributes and values. Attributes, are for each dataset
# as we saw above
#

# for all organisms

for (i in 1:length(org)){ # foreach organism
  print (org[i])
  for (j in 1:length(versions)){ # foreach of the 7 versions
  
    # load ensembl biomart in a variable named ensembl
  	#ensembl = useMart("ensembl", dataset=org[i]) # This is for the latest version of Ensembl - or you can get the website and use that as used below
  	ensembl = useMart(host=versions[j],biomart="ENSEMBL_MART_ENSEMBL", dataset=org[i]) # This is for the archived version (e.g. Ensembl v82 here)
  	
  	# generate file name
  	filename = paste(org[i], versions[j], sep = "_")
  	filename = paste("1_Paralogs_from_7_latest_versions", filename, sep="/")
  	filename = gsub(".org", ".txt", filename)
    
    # Get the attributes that I want to get for this organism
  	orgstring = strsplit(org[i], "_")[[1]] # get the name of organism
  	#print(orgstring[1]) # To test the name
   
    #attr = c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene", "hsapiens_paralog_subtype", "hsapiens_paralog_paralogy_confidence") # This is what I want to generate
    # generate the attribute array for the organism
    attr = c("ensembl_gene_id", 
             paste(orgstring[1], "_paralog_ensembl_gene", sep = ""),
             paste(orgstring[1], "_paralog_subtype", sep = ""),
             paste(orgstring[1], "_paralog_paralogy_confidence", sep = ""))
    
  	# var = getBM(attributes = attr, filters = 'biotype', values = 'protein_coding', mart = ensembl) # Earlier I was filtering PC genes, now I am not
  	var = getBM(attributes = attr, mart = ensembl)
    
  	write.table(var, file = filename, sep ="\t", quote=F, row.names = F)
  	
    rm("ensembl", "filename", "orgstring", "var", "attr") # remove some variables for the next run
    
  	print (versions[j])
  }	
}

