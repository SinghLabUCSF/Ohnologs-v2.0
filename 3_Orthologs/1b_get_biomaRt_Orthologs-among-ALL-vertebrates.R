# To get orthologs among all the 29 WGD organisms
####################################################################################


# load the biomart libraray
library("biomaRt")

# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read organism names for either 2R or 3R WGD. These have the Ensembl style names 
# that I need to get stuff from Biomart

#org = readLines("speciestest.txt")
#outgp = readLines("speciestest2.txt")

# For this version, the order doesn't matter.
org = readLines("all_vertebrate_species_ordered.txt") 

# define a filename tag
tag = "orthologs_Ensembl84.txt"

# load the biomart libraray
#library("biomaRt")



# for all 2R or 3R WGD organisms
for (i in 1:length(org)){
  
  # load ensembl biomart in a variable named ensembl for the current organism dataset=org[i]
  #ensembl = useMart("ensembl", dataset=org[i]) # This is for the latest version of Ensembl which is v84
  ensembl = useMart(host="mar2016.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", dataset=org[i])
  
  # get the name of WGD organism
  orgstring = strsplit(org[i], "_")[[1]]
  #print(orgstring[1]) # To test the name
  
  # For each of the outgroup
  for (j in 1:length(org)){
    
    # get the name of OUTGROUP organism
    ogpstr = strsplit(org[j], "_")[[1]]
    #print(ogpstr[1]) # To test the name
    
    # See if it's the same organism. If yes, skip and go to next
    if (ogpstr[1] == orgstring[1]){
      #print (paste ("***", ogpstr[1], " ", orgstring[1]))
      next;
    }
    
    # generate file name to be saved for this pair
    filename = paste(orgstring[1], ogpstr[1], tag, sep = "_")

    # Check if the file already exists, from previous run. If yes, skip and go to next. 
    # ************* UNCOMMENT FOR NEXT RUN ***************
    # if (file.exists(paste ("Vertebrate_Orthologs", "/", filename, sep = ""))){
    #   #print (filename)
    #   next;
    # }
    
    # If not print the name 
    print(filename)
    
    # Generate and attribute array for the organism - this would be the orthologs wrt the current outgroup organism
    attr = c("ensembl_gene_id",
             paste(ogpstr[1], "_homolog_ensembl_gene", sep = ""),
             paste(ogpstr[1], "_homolog_orthology_confidence", sep = ""))

    #print (attr) #test print

    # get the data from Biomart
    var = getBM(attributes = attr, mart = ensembl)

    # Save in a table
    write.table(var, file = paste ("Vertebrate_Orthologs", "/", filename, sep = ""), sep ="\t", quote=F, row.names = F)

    # print the pair that is just done
    print (c(org[i], org[j]))

    # remove some variables for the next run
    rm("ogpstr", "filename", "attr", "var")
    
  }
  # remove some variables for the next run
  rm("ensembl", "orgstring")
  
}

rm(list =ls())
