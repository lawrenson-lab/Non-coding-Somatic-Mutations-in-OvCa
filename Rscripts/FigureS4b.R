rm(list = ls())

library(cowplot)
library(ComplexHeatmap)
library(DESeq2)
library(DiffBind)
library(MASS)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(VennDiagram)
library(biomaRt)
library(circlize)
library(edgeR)
library(ggplot2)
library(org.Hs.eg.db)
library(plyr)
library(reshape2)
library(rtracklayer)
library(scales)

source('Rscripts/RegisterCores.R')
source('Rscripts/utils.R')
source('Rscripts/layer_IO.R')
source('Rscripts/layer_business.R')

col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')

x <- getRawCounts(filterGenes = F)
normalizedCounts <- getNormalizedCounts(x = x[,7:ncol(x)], gene_length = x$Length, scale = F)

plotCPM <- function(hgnc_symbol){
  .normalizedCounts <- normalizedCounts
  .bm <- bm[bm$hgnc_symbol == hgnc_symbol,]
  
  .normalizedCounts <- as.data.frame(.normalizedCounts)
  .normalizedCounts$ensembl_gene_id <- substr(rownames(.normalizedCounts),0,15)
  .normalizedCounts <- .normalizedCounts[.normalizedCounts$ensembl_gene_id %in% .bm$ensembl_gene_id,]
  .x <- merge(.bm, .normalizedCounts, by = 'ensembl_gene_id', sort = F)
  .x <- .x[,c(1,4,10:ncol(.x))]
  .x <- melt(.x)
  .x$variable <- as.character(.x$variable)
  .x$histotype <- substr(.x$variable,0,nchar(.x$variable)-2)
  .x$hgnc_symbol <- factor(x = .x$hgnc_symbol, levels = as.character(unique(.bm$hgnc_symbol)))
  colnames(.x)[4] <- 'CPM'
  .x$histotype <- factor(x = .x$histotype, levels =  unique(.x$histotype))
  p <- ggplot(.x, aes(histotype, CPM, col=histotype))
  p + geom_boxplot(outlier.color = NA, show.legend = F) + geom_jitter(width = 0.1, size = 0.5, show.legend = F) +
    scale_y_continuous(limits = c(0,max(.x$CPM))) + ggtitle(hgnc_symbol)
}



p <- plotCPM('LINC00963') + 
  ylab('Gene expression (CPM)') + 
  scale_color_brewer(palette = 'Set1') +
  theme_classic(base_size = 6) + 
  theme(axis.title.x = element_blank(), plot.margin = margin(0,0,0,0,'in'), axis.text = element_text(colour = 'black'))


w <- 3.125
pdf(file = 'figures/FigureS4b.pdf', width =w/2, height = w/2, useDingbats = F)
print(p)
dev.off()


