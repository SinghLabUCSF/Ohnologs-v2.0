# read organism names
org = readLines("species.txt")
#org = readLines("speciestest.txt")

# define a filename tag
tag = "84_GO_20160703.txt"

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


# Get the attributes to be listed in a variable
#attr = c("ensembl_gene_id", "strand", "chromosome_name", "external_gene_name", "start_position", "end_position", "band", "gene_biotype", "description", "go_id", "name_1006")
#attr = c("ensembl_gene_id", "strand", "chromosome_name", "external_gene_name", "start_position", "end_position", "band", "gene_biotype", "description")
attr = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006")

# for all organisms

for (i in 1:length(org)){
	
	# load ensembl biomart in a variable named ensembl
	# if the organism is sea urchin  - use ensembl metazoa biomart
	
	if (identical (org[i], "spurpuratus_gene_ensembl")){
		ensembl = useMart("metazoa_mart_21", dataset="spurpuratus_eg_gene") # The dataset name is also slightly different in this case - I am not using it here
	}
	# else use regular biomart and names
	else {
		ensembl = useMart("ensembl", dataset=org[i]) # This is for the latest version of Ensembl which is v84
	  #ensembl = useMart(host="sep2015.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", dataset=org[i]) # This is for the archived version (e.g. Ensembl v82 here)
	
	}
	
	# generate file name
	filename = paste(org[i], tag, sep = "")
	
	# var = getBM(attributes = attr, filters = 'biotype', values = 'protein_coding', mart = ensembl) # Earlier I was filtering PC genes, now I am not
	var = getBM(attributes = attr, mart = ensembl)
  
	write.table(var, file = filename, sep ="\t", quote=F, row.names = F)
	
  rm("ensembl", "filename", "var") # remove some variables for the next run
  
	print (org[i])
	
}
