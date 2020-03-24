rm(list = ls())

library(ggrepel)
library(reshape2)
library(plyr)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(cowplot)

source('Rscripts/RegisterCores.R')

cutoff <- 0.5
cols <- 3

dirs <- list.dirs(path = 'results/metascapeSEs', recursive = F)
data.long <- foreach(i=1:length(dirs), .combine = rbind)%do%{
  dir <- dirs[i]
  file <- paste0(dir, '/Enrichment_GO/_FINAL_GO.csv')
  .x <- read.csv(file = file, header = T)
  .x <- .x[match(x = unique(.x$BestEnrichmentInGroup), table = .x$BestEnrichmentInGroup),]
  return(.x)
}

data.long$X.GeneInHitList <- factor(data.long$X.GeneInHitList)
data.long$geneSet <- NA
data.long$geneSet[data.long$X.GeneInHitList == 86] <- 'MOC-specific'
data.long$geneSet[data.long$X.GeneInHitList == 99] <- 'CCOC-specific'
data.long$geneSet[data.long$X.GeneInHitList == 109] <- 'HGSOC-specific'
data.long$geneSet[data.long$X.GeneInHitList == 148] <- 'common'
data.long$label <- data.long$Description

data.long$LogP <- -data.long$LogP

plotPathways <- function(geneSet = 'CCOC-specific_pos', fill = '#008B8B'){
  data <- data.long[data.long$geneSet == geneSet,]
  data <- ddply(data, .(Description), summarize, LogP = max(LogP))
  data <- data[order(data$LogP, decreasing = T),]
  data$Description <- factor(data$Description, levels = data$Description)
  print(nrow(data))
  data <- data[1:(min(20,nrow(data))),]
  p <- ggplot(data = data, mapping = aes(x = Description, y = LogP)) + 
    geom_col(fill = fill) + 
    geom_hline(yintercept = -log10(0.05), linetype = 3) +
    theme_classic(base_size = 6) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}

p1 <- plotPathways(geneSet = 'CCOC-specific', fill = '#e41a1c')
p2 <- plotPathways(geneSet = 'HGSOC-specific', fill = '#4daf4a')
p3 <- plotPathways(geneSet = 'MOC-specific', fill = '#984ea3')
p4 <- plotPathways(geneSet = 'common', fill = '#999999')

pdf(file = 'figures/FigureS4d.pdf', width = 180/25.4, height = 180/25.4/2)
plot_grid(p1,p2,p3,p4, nrow = 1, ncol = 4, align = 'hv')
dev.off()
