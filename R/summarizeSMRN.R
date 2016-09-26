 # Make some summarizing plots for the SMRN with ggplot

rm(list = ls())
setwd('/Users/ereznik/Documents/small-molecule-regulation/R/')
library(ggplot2)
library(reshape2)
library(ggExtra)
library(ggrepel)

d = read.csv('../cache/iJO1366_SMRN.csv',header = TRUE,row.names = 1)

# Remove redundancy
d = d[,which(colnames(d) != 'bigg.reaction')]
d = d[!duplicated(d),]

act = d[d['Value']==1,]
inh = d[d['Value']==-1,]

uqrxns = unique(d$EC_number)
uqmets = unique(d$bigg.metabolite)
rxn = data.frame( 'Activating' = numeric(0), 'Inhibiting' = numeric(0),'Name' = character(0),
                  stringsAsFactors = FALSE)
met = rxn

##############################################################################
# Reactions
##############################################################################
for (rname in uqrxns){
  # Find number of activating/inhibiting edges
  rxn[rname,'Activating'] = length(which(act$EC_number == rname))
  rxn[rname,'Inhibiting'] = length(which(inh$EC_number == rname))
  rxn[rname,'Name'] = rname
}

# Remove labels for entities that are not outliers
rxn[rowSums(rxn[,c('Activating','Inhibiting')]) < 10,'Name'] = ''

rxnplot = ggplot(rxn,aes(Activating,Inhibiting,label = Name)) + stat_binhex(binwidth = c(2,0.5)) + 
  theme_bw(base_size = 18) + 
  scale_fill_gradient2(low = 'white',high = 'red',midpoint = 1,mid = 'gray') + 
  theme(legend.position="None") + 
  geom_text_repel(data = subset(rxn[rxn$Name != '',]),nudge_x = -3) + 
  xlab('Number of Activating Interactions') + ylab('Number of Inhibiting Interactions')
pdf(file = '../res/SMRN_Reaction_Histogram.pdf',onefile=FALSE)
ggMarginal(rxnplot,type = 'histogram')
dev.off()

##############################################################################
# Metabolites
##############################################################################
for (mname in uqmets){
  
  # Find number of activating/inhibiting edges
  met[mname,'Activating'] = length(which(act$bigg.metabolite == mname))
  met[mname,'Inhibiting'] = length(which(inh$bigg.metabolite == mname))
  met[mname,'Name'] = mname
}

# Remove labels
met[rowSums(met[,c('Activating','Inhibiting')]) < 7,'Name'] = ''

rxnplot = ggplot(met,aes(Activating,Inhibiting,label = Name)) + stat_binhex(binwidth = c(1,1)) + 
  theme_bw(base_size = 18) + 
  scale_fill_gradient2(low = 'white',high = 'red',midpoint = 1,mid = 'gray') + 
  theme(legend.position="None") + 
  geom_text_repel(data = subset(met[met$Name != '',]),nudge_x = -3) + 
  xlab('Number of Activating Interactions') + ylab('Number of Inhibiting Interactions')
pdf(file = '../res/SMRN_Metabolite_Histogram.pdf',onefile=FALSE)
ggMarginal(rxnplot,type = 'histogram')
dev.off()