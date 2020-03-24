rm(list = ls())

library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(scales)
library(UpSetR)

source('Rscripts/RegisterCores.R')

path <- 'results/RoseRes'
dirs <- list.dirs(path = path, recursive = F)
selectedDirs <- c(grep(pattern = 'ClearCell_241', x = dirs),
                   grep(pattern = 'ClearCell_511', x = dirs),
                   grep(pattern = 'ClearCell_3172', x = dirs),
                   grep(pattern = 'ClearCell_3547', x = dirs),
                   grep(pattern = 'ClearCell_3588', x = dirs),
                   grep(pattern = 'Endometrioid_291', x = dirs),
                   grep(pattern = 'Endometrioid_381', x = dirs),
                   grep(pattern = 'Endometrioid_589_Tu', x = dirs),
                   grep(pattern = 'Endometrioid_630', x = dirs),
                   grep(pattern = 'Endometrioid_703', x = dirs),
                   grep(pattern = 'HGSerous_229', x = dirs),
                   grep(pattern = 'HGSerous_270', x = dirs),
                   grep(pattern = 'HGSerous_429', x = dirs),
                   grep(pattern = 'HGSerous_550', x = dirs),
                   grep(pattern = 'HGSerous_561', x = dirs),
                   grep(pattern = 'Mucinous_230', x = dirs),
                   grep(pattern = 'Mucinous_380', x = dirs),
                   grep(pattern = 'Mucinous_464', x = dirs),
                   grep(pattern = 'Mucinous_652', x = dirs),
                   grep(pattern = 'Mucinous_3724', x = dirs))
title <- c(paste0('CCOC-',1:5),
           paste0('EnOC-',1:5),
           paste0('HGSOC-',1:5),
           paste0('MOC-',1:5))

SEAGs.union <- sort(unique(foreach(i=selectedDirs, .combine = c)%do%{
  title <- strsplit(x = dirs[i], split = '/')[[1]][3]
  print(title)
  x.SEs <- read.table(file = paste(dirs[i],list.files(path = dirs[i], pattern = "*_SuperEnhancers_ENHANCER_TO_GENE.ivyAnnotation.txt"), sep = '/'), header = T, sep = '\t')
  return(as.character(x.SEs$hgnc_symbol))
}))
SEAGs.ccoc <- sort(unique(foreach(i=selectedDirs[1:5], .combine = c)%do%{
  title <- strsplit(x = dirs[i], split = '/')[[1]][3]
  print(title)
  x.SEs <- read.table(file = paste(dirs[i],list.files(path = dirs[i], pattern = "*_SuperEnhancers_ENHANCER_TO_GENE.ivyAnnotation.txt"), sep = '/'), header = T, sep = '\t')
  return(as.character(x.SEs$hgnc_symbol))
}))
SEAGs.enoc <- sort(unique(foreach(i=selectedDirs[6:10], .combine = c)%do%{
  title <- strsplit(x = dirs[i], split = '/')[[1]][3]
  print(title)
  x.SEs <- read.table(file = paste(dirs[i],list.files(path = dirs[i], pattern = "*_SuperEnhancers_ENHANCER_TO_GENE.ivyAnnotation.txt"), sep = '/'), header = T, sep = '\t')
  return(as.character(x.SEs$hgnc_symbol))
}))
SEAGs.hgsoc <- sort(unique(foreach(i=selectedDirs[11:15], .combine = c)%do%{
  title <- strsplit(x = dirs[i], split = '/')[[1]][3]
  print(title)
  x.SEs <- read.table(file = paste(dirs[i],list.files(path = dirs[i], pattern = "*_SuperEnhancers_ENHANCER_TO_GENE.ivyAnnotation.txt"), sep = '/'), header = T, sep = '\t')
  return(as.character(x.SEs$hgnc_symbol))
}))
SEAGs.moc <- sort(unique(foreach(i=selectedDirs[16:20], .combine = c)%do%{
  title <- strsplit(x = dirs[i], split = '/')[[1]][3]
  print(title)
  x.SEs <- read.table(file = paste(dirs[i],list.files(path = dirs[i], pattern = "*_SuperEnhancers_ENHANCER_TO_GENE.ivyAnnotation.txt"), sep = '/'), header = T, sep = '\t')
  return(as.character(x.SEs$hgnc_symbol))
}))

data <- as.data.frame(ifelse(test = cbind(SEAGs.union %in% SEAGs.ccoc,
                                          SEAGs.union %in% SEAGs.enoc,
                                          SEAGs.union %in% SEAGs.hgsoc,
                                          SEAGs.union %in% SEAGs.moc), yes = 1, no = 0))
colnames(data) <- c('CCOC', 'EnOC', 'HGSOC', 'MOC')
rownames(data) <- SEAGs.union

plot.size <- 188/25.4
colors <- brewer.pal(n = 4, name = 'Set1')
pdf(file = 'figures/Figure2c.pdf', width = plot.size/2, height = plot.size/2, pointsize = 6, onefile = F)
upset(data = data, 
      order.by = 'freq', 
      mainbar.y.label = 'Number of genes', 
      sets.bar.color = colors, 
      keep.order = T, 
      sets = c('CCOC', 'EnOC', 'HGSOC', 'MOC'), 
      sets.x.label = 'Number of genes',
      main.bar.color = c('black', colors[1], colors[3], colors[4], 'black', colors[2], rep('black',9)))
dev.off()

