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

dirs <- list.dirs(path = 'results/metascape', recursive = F)
data.long <- foreach(i=1:length(dirs), .combine = rbind)%do%{
  dir <- dirs[i]
  file <- paste0(dir, '/Enrichment_GO/_FINAL_GO.csv')
  .x <- read.csv(file = file, header = T)
  .x <- .x[match(x = unique(.x$BestEnrichmentInGroup), table = .x$BestEnrichmentInGroup),]
  colnames(.x)[1:2] <- c("X_MEMBER_Input.ID", "X_LogP_Input.ID"  )
  return(.x)
}

data.long$X.GeneInHitList <- factor(data.long$X.GeneInHitList)
data.long$geneSet <- NA
data.long$geneSet[data.long$X.GeneInHitList == 100] <- 'MOC-specific_pos'
data.long$geneSet[data.long$X.GeneInHitList == 105] <- 'HGSOC-specific_pos'
data.long$geneSet[data.long$X.GeneInHitList == 121] <- 'CCOC-specific_pos'
data.long$geneSet[data.long$X.GeneInHitList == 182] <- 'MOC-specific_neg'
data.long$geneSet[data.long$X.GeneInHitList == 276] <- 'CCOC-specific_neg'
data.long$geneSet[data.long$X.GeneInHitList == 620] <- 'HGSOC-specific_neg'
data.long$geneSet[data.long$X.GeneInHitList == 2979] <- 'EOC-common'
data.long$label <- data.long$Description
data.long$label[data.long$LogP > log10(0.05) | data.long$Enrichment < 50] <- NA

data.long$LogP <- -data.long$LogP

plotPathways <- function(geneSet = 'CCOC-specific_pos', fill = '#008B8B'){
  data <- data.long[data.long$geneSet == geneSet,]
  data <- ddply(data, .(Description), summarize, LogP = max(LogP))
  data <- data[order(data$LogP, decreasing = T),]
  data$Description <- factor(data$Description, levels = data$Description)
  print(nrow(data))
  data <- data[1:(min(21,nrow(data))),]
  p <- ggplot(data = data, mapping = aes(x = Description, y = LogP)) + 
    geom_col(fill = fill) + 
    geom_hline(yintercept = -log10(0.05), linetype = 3) +
    theme(text = element_text(family = 'sans', size = 5), 
          panel.grid = element_blank(), panel.background = element_blank(), 
          axis.text = element_text(family = 'sans', size = 5),
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0, 1, 0, 8), "mm"))
  return(p)
}

p1 <- plotPathways(geneSet = 'CCOC-specific_pos', fill = '#e41a1c')
p2 <- plotPathways(geneSet = 'CCOC-specific_neg', fill = '#fbb4ae')
p5 <- plotPathways(geneSet = 'HGSOC-specific_pos', fill = '#4daf4a')
p6 <- plotPathways(geneSet = 'HGSOC-specific_neg', fill = '#ccebc5')
p7 <- plotPathways(geneSet = 'MOC-specific_pos', fill = '#984ea3')
p8 <- plotPathways(geneSet = 'MOC-specific_neg', fill = '#decbe4')
p9 <- plotPathways(geneSet = 'EOC-common', fill = '#999999')

pdf(file = 'figures/Figure1k.pdf', width = 188/25.4, height = 188/25.4/2)
plot_grid(p1,p5,p7, nrow = 1, ncol = 3, align = 'hv')
dev.off()

pdf(file = 'figures/FigureS3.pdf', width = 2*3.125, height = 3*3.125)
plot_grid(p6,p9,p2, p8, nrow = 2, ncol = 2, align = 'hv')
dev.off()




