# To get orthologs between WGD organisms and their outgroups
####################################################################################
#
# 2R-WGD organisms: 29
# 2R WGD outgroups in Ensembl: 4 (amphioxus is separate)
#
# 3R-WGD organisms: 4 fish species
# 3R WGD outgroups in Ensembl: 6
####################################################################################

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read organism names for either 2R or 3R WGD. These have the Ensembl style names 
# that I need to get stuff from Biomart

# ---------------- for test -------------------
#org = readLines("speciestest.txt")
#outgp = readLines("speciestest2.txt")

# ---------------- for 2R WGD -------------------
#org = readLines("2r_wgd_species.txt")
#outgp = readLines("2r_outgroup_species.txt")

# ---------------- for 3R WGD -------------------
org = readLines("3r_wgd_species.txt")
outgp = readLines("3r_outgroup_species.txt")

# define a filename tag
tag = "orthologs_Ensembl84.txt"

# load the biomart libraray
library("biomaRt")


# for all 2R or 3R WGD organisms
for (i in 1:length(org)){
	
	# load ensembl biomart in a variable named ensembl for the current organism dataset=org[i]
  #ensembl = useMart("ensembl", dataset=org[i]) # This is for the latest version of Ensembl which is v84
  ensembl = useMart(host="mar2016.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", dataset=org[i])
  
	# get the name of WGD organism
	orgstring = strsplit(org[i], "_")[[1]]
	#print(orgstring[1]) # To test the name
	
  # For each of the outgroup
	for (j in 1:length(outgp)){
  
	  # get the name of OUTGROUP organism
	  ogpstr = strsplit(outgp[j], "_")[[1]]
	  #print(ogpstr[1]) # To test the name
    
    # generate file name to be saved for this pair
	  filename = paste(orgstring[1], ogpstr[1], tag, sep = "_")

    # Generate and attribute array for the organism - this would be the orthologs wrt the current outgroup organism
    attr = c("ensembl_gene_id", 
             paste(ogpstr[1], "_homolog_ensembl_gene", sep = ""),
             paste(ogpstr[1], "_homolog_orthology_confidence", sep = ""))
	
    #print (attr) #test print
    
    # get the data from Biomart
    var = getBM(attributes = attr, mart = ensembl)
    
    # Save in a table 
    write.table(var, file = filename, sep ="\t", quote=F, row.names = F)
	
    # print the pair that is just done
  	print (c(org[i], outgp[j]))
    
	  # remove some variables for the next run
	  rm("ogpstr", "filename", "attr", "var")
    
	}
  # remove some variables for the next run
  rm("ensembl", "orgstring")
  
}

rm(list =ls())



