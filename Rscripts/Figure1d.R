rm(list = ls())

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

#chip-seq
gr <- getAllHistotypeSpecificH3K27acSites()
gr.consensus <- getH3K27acTumorConsensus()

#rna-seq
x <- getRawCounts()
deGenes <- getDEGenes(x = x)
normalizedCounts <- getNormalizedCounts(x[,7:ncol(x)],T)

w <- 88/25.4
pdf(file = paste0('figures/Figure1d.pdf'), width = w, height = w, onefile = T, useDingbats = F)
dMax <- 1000000
histotypes <- c('CCOC','EnOC','HGSOC','MOC')
for(fold in c('pos','neg')){
  for(histotype.1 in histotypes){
    print(histotype.1)
    df <- foreach(histotype.2=histotypes, .combine = rbind) %dopar% {
      cpm <- foreach(th=seq(-dMax,dMax,1000), .combine = c) %dopar% {
        if(fold == 'pos'){
          .gr <- gr[gr$histotype == histotype.1 & gr$Fold>0]  
        }else{
          .gr <- gr[gr$histotype == histotype.1 & gr$Fold<0]
        }
        if(th<0){
          start(.gr) <- start(.gr)+th
        }else if(th>0){
          end(.gr) <- end(.gr)+th
        }
        .cpm <- mean(normalizedCounts[substr(rownames(normalizedCounts),0,15) %in% bm$ensembl_gene_id[bm$entrezgene %in% subsetByOverlaps(genes,.gr)$gene_id],startsWith(colnames(normalizedCounts),histotype.2)])
        return(.cpm)
      }
      .df <- data.frame(histotype.1 = rep(histotype.1, length(cpm)),
                        histotype.2 = rep(histotype.2, length(cpm)),
                        th = seq(-dMax,dMax,1000)/1000,
                        cpm = cpm)
      return(.df)
    }
    ylim <- max(abs(df$cpm))
    print(ylim)
    ylim <- 0.72
    p <- ggplot(data = df, mapping = aes(th,cpm, col=histotype.2, linetype=histotype.1))
    print(
      p + geom_line(show.legend = F, size = 2) + scale_color_manual(values = col) + scale_x_continuous(labels = comma) + 
        scale_y_continuous(limits = c(-ylim,ylim)) + theme_light() +
        theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
        geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
    )
  }
}
dev.off()