# Read in species names and extract taxonomic information.
# The queries sometimes ask for user input, be prepared to sit and wait. May take a few hours...:(
# We built in some functionality to "restart" if a _temp file exists and doRestart = TRUE

setwd('/Users/ereznik/Documents/small-molecule-regulation/R/')
library(taxize)
library(data.table)

doRestart = TRUE #if TRUE, then we read in the _temp file we had been working on before

simpleCap <- function(x) {
  paste(toupper(substring(x, 1,1)), substring(x, 2),
        sep="", collapse=" ")
}

# Read in taxonomic information for each file in BRENDA query
brendadir = '../data/brenda_query_2016_06_08/'
species = c()
for (f in list.files(brendadir)){
  temp = fread(paste(brendadir,f,sep = ''),header = TRUE,data.table = FALSE,sep = ',',encoding="Latin-1")
  species = c(species,temp[,2])
}
uqspecies = unique(species)
#uqspecies = sapply(uqspecies,simpleCap)

# Resolve species names
starttime = proc.time()

if (doRestart){
  res = read.csv('../cache/TaxonomicData_temp.csv',header = TRUE,row.names = 1)
}else{
  res = data.frame(matrix(0,0,8))
  colnames(res) = c('kingdom','phylum','class','order','family','genus','species','RevisedName')
}


for (i in 1:length(uqspecies)){
  specname = uqspecies[i]
  
  if (doRestart & specname %in% rownames(res)){
    print(paste('Skipping',specname,'already imported...'))
    next # we already imported this one
  }
  
  # Initialize the rows with NA
  res[specname,] = NA

  # We may have problems modifying the species name, we can fix this
  newspecname = tryCatch({
    species_resolve = gnr_resolve( tolower(specname) )
    newspecname = species_resolve[1,'matched_name'] # Dumbly select the first row among all the options
  }, error = function(e) {
    print(paste('Skipping the renaming of species:',specname))
    newspecname = NULL
  },warning = function(war) {
    print(paste('Skipping the renaming of species:',specname))
    newspecname = NULL
  }) #end tryCatch
  
  # If the new species name is null, that means we had an issue, skip it
  if (is.null(newspecname)){
    next
  }
  
  # Get the taxonomic information if we've made it this far
  species_classify = classification(newspecname, db = 'itis')
  
  temp = species_classify[[1]]
  if( is.null(dim(temp)) ) {
    print(paste('Skipping',specname,'unable to classify.'))
    next
  }else{
    rownames(temp) = temp$rank
    
    writedata_types = intersect( rownames(temp), colnames(res) )
    writedata = temp[ writedata_types,'name' ]
    res[specname,writedata_types] = writedata
    res[specname,'RevisedName'] = newspecname
  }
  if (i%%10 == 0){
    print(i)
    # Save the temporary data in case of a crash
    write.csv(res,'../cache/TaxonomicData_temp.csv')
  }

}

print(proc.time() - starttime)
write.csv(res,'../cache/TaxonomicData.csv')
