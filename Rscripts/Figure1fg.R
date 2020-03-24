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

w <- 89/2.54/10
col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')

gr <- getAllHistotypeSpecificH3K27acSites()
gr.consensus <- getH3K27acTumorConsensus()
x <- getRawCounts()
deGenes <- getDEGenes(x)
deGenes.pos <- getDEGenes(x = x, positiveFold = T, negativeFold = F)
deGenes.neg <- getDEGenes(x = x, positiveFold = F, negativeFold = T)
normalizedCounts <- getNormalizedCounts(x[,7:ncol(x)], x$Length, F)

external_gene_name <- 'WFDC2'

source('Rscripts/mappingUtils.R')
library(gridExtra)
mapping.init()
df <- mapEnhancerToGene(gr.enhancer = GRanges(seqnames = 'chr20', ranges = IRanges(start = 44095981, end = 44101008)), 
                        ensembl_gene_id = 'ENSG00000101443', type = 'equal')
cor.test(x = df$x, y = df$y, method = 'spearman', alternative = 'greater')

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
    scale_y_continuous(limits = c(0,max(.x$CPM)))
}


p1 <- plotCPM(external_gene_name)
p1 <- p1 + 
  ylab('WFDC2 expression\nRNA-seq (CPM)') + 
  scale_color_brewer(palette = 'Set1') +
  theme_classic(base_size = 6) + 
  theme(axis.title.x = element_blank(), plot.margin = margin(0,0,0,0,'in'), axis.text = element_text(colour = 'black'))


p2 <- ggplot(data = df, mapping = aes(x, y)) + 
  geom_point(show.legend = F, mapping = aes(col=histotype), size = 0.5) + 
  geom_smooth(method=lm, se = F, col = 'black', size = 0.5) + 
  annotate(geom = 'text', x = 512, y = c(2,1), label = c('r = 0.78',expression(paste('p = 5.5 x 10'^'-5'))), size = 2, hjust = 0, vjust = 0) + 
  scale_color_brewer(palette = 'Set1') + 
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2') + 
  xlab('WFDC2 promoter activity\nH3K27ac ChIP-seq\n(normalized read count)') + 
  ylab('WFDC2 expression\nRNA-seq (CPM)') + 
  theme_classic(base_size = 6) +
  theme(plot.margin = margin(0,0,0,0,'in'), axis.text = element_text(colour = 'black'))


w <- 3.125
pdf(file = 'figures/Figure1fg.pdf', width = w/2, height = w, useDingbats = F)
plot_grid(p1,p2,align = 'hv', nrow = 2, ncol = 1, labels = c('f','g'), label_size = 10)
dev.off()


