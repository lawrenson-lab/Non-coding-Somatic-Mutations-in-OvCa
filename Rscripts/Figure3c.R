library(GenomicRanges)
rm(list = ls())
x1 <- read.table(file = 'results/aim1/aim1_FMREs_CCOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x2 <- read.table(file = 'results/aim1/aim1_FMREs_EnOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x3 <- read.table(file = 'results/aim1/aim1_FMREs_HGSOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x4 <- read.table(file = 'results/aim1/aim1_FMREs_CCOC+EnOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x5 <- read.table(file = 'results/aim1/aim1_FMREs_CCOC+EnOC+HGSOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x1 <- GRanges(x1)
x2 <- GRanges(x2)
x3 <- GRanges(x3)
x4 <- GRanges(x4)
x5 <- GRanges(x5)
x1 <- x1[x1$nSamples > 2 & x1$p.adjust <= 0.25]
x2 <- x2[x2$nSamples > 2 & x2$p.adjust <= 0.25]
x3 <- x3[x3$nSamples > 2 & x3$p.adjust <= 0.25]
x4 <- x4[x4$nSamples > 2 & x4$p.adjust <= 0.25]
x5 <- x5[x5$nSamples > 2 & x5$p.adjust <= 0.25]
x1 <- x1[order(x1$p.value)]
x2 <- x2[order(x2$p.value)]
x3 <- x3[order(x3$p.value)]
x4 <- x4[order(x4$p.value)]
x5 <- x5[order(x5$p.value)]
x1$histotype <- 'CCOC'
x2$histotype <- 'EnOC'
x3$histotype <- 'HGSOC'
x4$histotype <- 'CCOC+EnOC'
x5$histotype <- 'ALL'
x <- c(x1,x2,x3,x4,x5)
x <- reduce(x)

x1 <- read.table(file = 'results/aim1/aim1_FMREs_CCOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x2 <- read.table(file = 'results/aim1/aim1_FMREs_EnOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x3 <- read.table(file = 'results/aim1/aim1_FMREs_HGSOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x4 <- read.table(file = 'results/aim1/aim1_FMREs_CCOC+EnOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x5 <- read.table(file = 'results/aim1/aim1_FMREs_CCOC+EnOC+HGSOC_minOverlap4_flanklin250bp.txt', header = T, sep = '\t', quote = '')
x1 <- GRanges(x1)
x2 <- GRanges(x2)
x3 <- GRanges(x3)
x4 <- GRanges(x4)
x5 <- GRanges(x5)

source('Rscripts/RegisterCores.R')
mcols(x) <- foreach(i=1:length(x), .combine = rbind)%dopar%{
  c(ifelse(test = length(subsetByOverlaps(x1, x[i])) == 0, no = min(subsetByOverlaps(x1, x[i])$p.value), yes = NA),
    ifelse(test = length(subsetByOverlaps(x2, x[i])) == 0, no = min(subsetByOverlaps(x2, x[i])$p.value), yes = NA),
    ifelse(test = length(subsetByOverlaps(x3, x[i])) == 0, no = min(subsetByOverlaps(x3, x[i])$p.value), yes = NA),
    ifelse(test = length(subsetByOverlaps(x4, x[i])) == 0, no = min(subsetByOverlaps(x4, x[i])$p.value), yes = NA),
    ifelse(test = length(subsetByOverlaps(x5, x[i])) == 0, no = min(subsetByOverlaps(x5, x[i])$p.value), yes = NA))
}
x <- as.data.frame(x)
colnames(x)[6:10] <- c('CCOC', 'EnOC', 'HGSOC', 'CCOC+EnOC', 'ALL')
library(reshape2)
x$RE <- paste0(x$seqnames, ':', x$start, '-', x$end)
x$RE[x$RE == 'chr3:126422837-126423886'] <- 'CHCHD6 promoter'
x$RE[x$RE == 'chr5:168727067-168728477'] <- 'SLIT3 promoter'
x$RE[x$RE == 'chr10:3827366-3829505'] <- 'KLF6 promoter'
x$RE[x$RE == 'chr16:22307741-22310284'] <- 'POLR3E promoter'
x$RE[x$RE == 'chr16:46863901-46865831'] <- 'C16orf87 promoter'
x$RE[x$RE == 'chr19:40747862-40749385'] <- 'AKT2 promoter'
x$RE[x$RE == 'chrX:15692029-15694894'] <- 'CA5BP1 promoter'
x$RE[x$RE == 'chrX:74144064-74145465'] <- 'NEXMIF promoter'
x$RE[x$RE == 'chrX:153990093-153992603'] <- 'DKC1 promoter'
x$RE[x$RE == 'chrX:53709258-53711891'] <- 'HUWE1 promoter'

library(ComplexHeatmap)
library(circlize)
matrix <- x[,c(6,7,8)]
matrix[matrix>0.05] <- 1
rownames(matrix) <- x$RE
matrix <- -log10(matrix)
matrix[is.na(matrix)] <- 0
matrix <- matrix[rowSums(matrix)>0,]
w <- 180/25.4
pdf(file = 'figures/Figure3c.pdf', width = 4*w/10, height = 4*w/10, onefile = T)
Heatmap(matrix = matrix, col = colorRamp2(c(0, max(matrix, na.rm = T)), c("#ffffff", '#008b8b')), 
        name = '-log10(p.value)', cluster_columns = F, row_names_max_width = unit(5, 'in'), na_col = "#ffffff", show_row_dend = T,
        cluster_rows = T, clustering_distance_rows = 'canberra', clustering_method_rows = 'ward.D2', 
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 6, fontface = "bold")))
dev.off()
