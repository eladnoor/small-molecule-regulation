# Visualize cross-species interactions

rm(list = ls())
setwd('/Users/ereznik/Documents/small-molecule-regulation/R/')
library(ggplot2)
library(pheatmap)
library(data.table)

# Issues:
# 1. In the original python script, we are removing rows and columns with sum less than 3. We are repeating here. 
# It should make only a small difference.

# Function for cosine similarity
cos.sim = function(ix,X) 
{
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}   

# Make plotting more restrictive if you like
minval = 70

# Drop certain EC numbers
dropEC = c('2.7.7.27')
dropligand = c('2')

inh = read.csv('../cache/inh_crosstab.csv',header = TRUE,row.names = 1,check.names = FALSE)
act = read.csv('../cache/act_crosstab.csv',header = TRUE,row.names = 1,check.names = FALSE)

if (length(dropEC)!=0){
  act = act[-which(rownames(act) %in% dropEC),]
  inh = inh[-which(rownames(inh) %in% dropEC),]
}

if (length(dropligand)!=0){
  act = act[,-which(colnames(act) %in% dropligand)]
  inh = inh[,-which(colnames(inh) %in% dropligand)]
}

acttrim = act[-which(rowSums(act) < minval),-which(colSums(act)<minval)]
inhtrim = inh[-which(rowSums(inh) < minval),-which(colSums(inh)<minval)]

acttrim = as.matrix(acttrim)
inhtrim = as.matrix(inhtrim)
#pheatmap(act)
#pheatmap(inhtrim)

# Compute cosine similarity on full data
d2calc = inhtrim
n = dim(d2calc)[1]
cmb = expand.grid(i=1:n, j=1:n) 
cossim = matrix( apply(cmb,1,function(x){cos.sim(x,d2calc)}), n,n )
rownames(cossim) = rownames(d2calc)
pheatmap(cossim,filename = '/Users/ereznik/Documents/temp/tempheatmap.pdf',fontsize = 5)
