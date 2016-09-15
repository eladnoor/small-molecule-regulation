# Compare revised taxonomy code to old taxonomy code
setwd('/Users/ereznik/Documents/small-molecule-regulation/R/checks/')
d1 = read.csv('../../cache/TaxonomicData.csv',header = TRUE,row.names = 1)
d2 = read.csv('TaxonomicData_v1.csv',header = TRUE,row.names = 1)

# Find intersecting rownames
ixnames = intersect(rownames(d1),rownames(d2))
offdiag = c()
# For each of the intersecting names, make a table reporting the agreement of the annotation
for (col in c('kingdom','phylum','class')){
  temp = data.frame('New' = d1[ixnames,col], 'Old' = d2[ixnames,col])
  print(table(temp))
  offdiag[col] = sum(table(temp)) - sum(diag(table(temp)))
  
}