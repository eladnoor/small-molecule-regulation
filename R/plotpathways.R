# Summarize pathway data

rm(list = ls())
setwd('/Users/ereznik/Documents/small-molecule-regulation/R')
library(ggExtra)
library(ggplot2)
library(reshape2)
library(NMF)
library(pheatmap)

act = read.csv('../res/pathway_met_histograms_activating.csv',header = TRUE,row.names = 1,check.names = FALSE)
inh = read.csv('../res/pathway_met_histograms_inhibiting.csv',header = TRUE,row.names = 1,check.names = FALSE)
act[is.na(act)] = 0
inh[is.na(inh)] = 0
# Cluster the data
# for (dname in c('act','inh')){
#   d = get(dname)
#   
#   distmets = dist(as.matrix(d), method = "euclidean")
#   fitmets = hclust(distmets,method = 'ward.D')
#   clustmets = cutree(fitmets,k = 10)
#   ordmets = names(clustmets)[order(clustmets)]
#   
#   distpaths = dist(t(d),method = 'euclidean')
#   fitpaths = hclust(distpaths,method = 'ward.D')
#   clustpaths = cutree(fitpaths,k = 10)
#   ordpaths = names(clustpaths)[order(clustpaths)]
#   
#   d$Metabolite = rownames(d)
#   d2 = melt(d,id.vars = 'Metabolite')
#   d2$Metabolite = factor(d2$Metabolite,levels = ordmets)
#   d2$variable = factor(d2$variable,levels = ordpaths)
#   
#   
#   p1 = ggplot(d2,aes(Metabolite,variable,fill = value)) + geom_tile() + 
#     scale_fill_gradient2()
#   
# }

# Make colors
paletteLength = 50
myColor = colorRampPalette(c("white", "blue","red"))(paletteLength)

# Cluster with aheatmap, plot with pheatmap
inh2 = inh[which(rowSums(inh) > 7),]
act2 = act[which(rowSums(act) > 2),]
inhcluster = aheatmap(inh2)
actcluster = aheatmap(act2)
inhbreaks = c(seq(max(inh2)/paletteLength, max(inh2), length.out=floor(paletteLength)))
actbreaks = c(seq(max(act2)/paletteLength, max(act2), length.out=floor(paletteLength)))

pheatmap(inh2[inhcluster[[2]],inhcluster[[4]]],
         breaks = inhbreaks,color = myColor,
         cluster_rows = FALSE,cluster_cols = FALSE,
         filename = '../res/pathway_inhibiting.pdf',height = 8,width = 8)

pheatmap(act2[actcluster[[2]],actcluster[[4]]],
         breaks = actbreaks,color = myColor,
         cluster_rows = FALSE,cluster_cols = FALSE,
         filename = '../res/pathway_activating.pdf',height = 8,width = 8)

