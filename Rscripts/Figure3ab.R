rm(list = ls())

library(qqman)
library(GenomicRanges)

# parameters
histotype <- 'HGSOC'
minOverlap <- 4
flanking <- 250

gr.ROI <- GRanges(read.table(file = paste0('results/aim1/aim1_FMREs_',histotype,'_minOverlap',minOverlap,'_flanklin',flanking,'bp.txt'), header = T, sep = '\t', quote = ''))

source('Rscripts/layer_presentation.R')

w <- 2*4
pdf(file = paste0('figures/Figure3b.pdf'), width = w, height = w/2, useDingbats = F)
fdr <- 0.25
.x1 <- -log10(min(gr.ROI$p.value[gr.ROI$p.adjust > fdr]))/2
.x2 <- -log10(max(gr.ROI$p.value[gr.ROI$p.adjust <= fdr]))/2
print(plotManhattan(gr.H3K27ac = gr.ROI, pvalues = gr.ROI$p.value, suggestiveline = 100, genomewideline = 100))
dev.off()

w <- 4
pdf(file = paste0('figures/Figure3a.pdf'), width = 0.85*w, height = w, useDingbats = F)
nPass <- sum(gr.ROI$p.adjust <= fdr)
print(qq(pvector = sort(gr.ROI$p.value), col = c(rep('#4daf4a',nPass),rep('#ccebc5',length(gr.ROI$p.value)-nPass))))
dev.off()
