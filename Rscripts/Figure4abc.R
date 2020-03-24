rm(list = ls())

library(GenomicRanges)
library(ComplexHeatmap)
require(circlize)
library(reshape2)
library(ggrepel)
library(qqman)

source('Rscripts/RegisterCores.R')
source('Rscripts/layer_presentation.R')

histotypes <- c('HGSOC', 'CCOC', 'EnOC')
minOverlap <- 4
flanking <- 0

selectedGenes <- c('PAX8', 'ZNF217', 'ARID1A', 'HOXD9', 'MECOM', 'KLF6',
                   'MEIS1', 'SEPT9', 'TELO2', 'ZSCAN16', 'SOX17', 'WT1', 'PBX1', 'MUC16', 
                   'LAMA5', 'PTK2', 'LAMC2', 'WNT7A', 'CTBP1', 'CCDC6', 'PAX8',
                   'PBX1', 'PAX8', 'TMPRSS2', 'PTK2', 'ASPSCR1', 'MYC',
                   'MUC16', 'CRABP1', 'CRABP2', 'EPCAM', 'ER', 
                   'SPON1', 'SPON2', 'WFDC2', 'HNF1A', 'HNF1B', 
                   'IGF1', 'IGF2', 'CDH6', 'CDH10', 'MKI67', 
                   'KISS1', 'ST14', 'TMPRSS6', 'MSLN', 'MSLNL', 
                   'MIF', 'MMP7', 'MMP20', 'CDKN1A', 'TP53', 
                   'PAX8', 'PGR', 'SLPI', 'TACSTD2', 'WT1',
                   'BHLHE40', 'BHLHE41', 'BCL2L1', 'EMX2', 'ARID1B')
x.seags <- read.table(file = 'data/SEAG_union.txt', header = T, sep = '\t', quote = '')
min <- 3
selectedGenes <- x.seags$hgnc_symbol[x.seags$n.CCOC >= min & x.seags$n.EnOC >= min & x.seags$n.HGSOC >= min]


histotype <- 'HGSOC'
title <- paste0('aim2.EOC_',histotype, '.minOverlap', minOverlap)
x <- read.table(file = paste0('results/',title,'.',flanking,'bp.filt.txt'), header = T, sep = '\t', quote = "")
colnames(x)[6:7] <- c('start', 'end')
x$strand <- ifelse(x$strand == 1, '+', '-')
x$chromosome_name <- paste0('chr', x$chromosome_name)
gr <- GRanges(x)

w <- 2*4
pdf(file = 'figures/Figure4c.pdf', width = w, height = w/2, useDingbats = F)
fdr <- 0.25
.x1 <- -log10(min(gr$p.value[gr$p.adjust > fdr]))/2
.x2 <- -log10(max(gr$p.value[gr$p.adjust <= fdr]))/2
plotManhattan(gr.H3K27ac = gr, pvalues = gr$p.value)
dev.off()

w <- 4
pdf(file = 'figures/Figure4a.pdf', width = w, height = w, useDingbats = F)
nPass <- sum(gr$p.adjust <= 0.25)
print(qq(pvector = gr$p.value, col = c(rep('#4daf4a',nPass),rep('#ccebc5',length(gr$p.value)-nPass))))
dev.off()






hgnc_symbols <- foreach(i=1:length(histotypes), .combine = c)%dopar%{
  histotype <- histotypes[i]
  title <- paste0('aim2.EOC_',histotype, '.minOverlap', minOverlap)
  x <- read.table(file = paste0('results/',title,'.',flanking,'bp.filt.txt'), header = T, sep = '\t', quote = "")
  x$FC <- x$nSamples/x$exp
  x$significance <- ifelse(x$p.value <= 0.001, 'p.value < 0.001', 
                           ifelse(x$p.value <= 0.01, 'p.value < 0.01', 
                                  ifelse(x$p.value <= 0.05, 'p.value < 0.05', 'others')))
  x$significance <- factor(x$significance, levels = c('p.value < 0.001', 'p.value < 0.01', 'p.value < 0.05', 'others'))
  x <- x[x$nSamples >= 0.05*x$nTotSamples,]
  return(as.character(x$hgnc_symbol[rank(x$p.value) <= 10]))
}

hgnc_symbols <- sort(unique(c(hgnc_symbols, as.character(selectedGenes))))
hgnc_symbols <- hgnc_symbols[nchar(hgnc_symbols) > 0]
p.values <- foreach(i=1:length(histotypes), .combine = cbind)%dopar%{
  histotype <- histotypes[i]
  title <- paste0('aim2.EOC_',histotype, '.minOverlap', minOverlap)
  x <- read.table(file = paste0('results/',title,'.',flanking,'bp.filt.txt'), header = T, sep = '\t', quote = "")
  x <- x[x$nSamples >= 0.05*x$nTotSamples,]
  return(x$p.adjust[match(x = hgnc_symbols, table = x$hgnc_symbol)])
}

matrix <- -log10(p.values)
matrix[is.na(matrix) | matrix <= -log10(0.25)] <- 0
colnames(matrix) <- histotypes
rownames(matrix) <- hgnc_symbols

matrix <- matrix[,c(2,3,1)]
matrix <- matrix[rowSums(matrix)>0,]
pdf(file = 'figures/Figure4b.pdf', width = 2*180/25.4/7, height = 6*225/25.4/10, useDingbats = T)
Heatmap(matrix = matrix, col = colorRamp2(c(0, max(matrix, na.rm = T)), c("#ffffff", '#008b8b')), 
        name = '-log10(FRD)', cluster_columns = F, row_names_max_width = unit(5, 'in'), na_col = "#ffffff", show_row_dend = T, cluster_rows = T,
        row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7), row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
        clustering_distance_rows = 'canberra', clustering_method_rows = 'ward.D2')
dev.off()
