rm(list = ls())
library(poibin)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(biomaRt)
library(qqman)

source('Rscripts/RegisterCores.R')
source('Rscripts/creagUtils.R')
source('Rscripts/utils.R')
source('Rscripts/ConsensusREs.R')
source('Rscripts/layer_presentation.R')

getOCMapping <- function(){
  oc.seqnames <- paste0('chr',c(1:22,'X'))
  creag <- getCREAGs()
  mapping <- data.frame(seqnames = seqnames(creag), start = start(creag), end = end(creag), ensembl_gene_id = creag$ensembl_gene_id)
  gr.mapping <- GRanges(mapping)
  gr.promoters <- getPromoters(.ensembl)
  mcols(gr.promoters) <- data.frame(ensembl_gene_id = gr.promoters$ensembl_gene_id)
  gr.mapping <- c(gr.mapping, gr.promoters)
  gr.mapping <- gr.mapping[seqnames(gr.mapping) %in% oc.seqnames]
  seqlevels(gr.mapping) <- oc.seqnames
  gr.mapping <- sort(gr.mapping)
  return(gr.mapping)
}

load(file = 'RData/noncodingSNVs.RData')
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/wgEncodeDacMapabilityConsensusExcludable.bed')) == 0]
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/seq.cov1.ONHG19.bed.gz')) == 0]
all.SNVs <- gr.SNV


minOverlap <- 4
flanking <- 0
histotype <- 'HGSOC'
title <- paste0('aim2.EOC_',histotype, '.minOverlap', minOverlap)
creag.init(minOverlap = minOverlap, histotype = histotype)
gr.ROI <- getOCMapping()
ensembl_gene_ids <- sort(unique((split(gr.ROI, strand(gr.ROI))$'*')$ensembl_gene_id))
gr.ROI <- gr.ROI[gr.ROI$ensembl_gene_id %in% ensembl_gene_ids]

if(histotype == 'CCOC'){
  gr.SNV <- all.SNVs[endsWith(x = all.SNVs$dataset, suffix = 'CCOC')]
}else if(histotype == 'EnOC'){
  gr.SNV <- all.SNVs[endsWith(x = all.SNVs$dataset, suffix = 'ENOC')]
}else if(histotype == 'HGSOC'){
  gr.SNV <- all.SNVs[endsWith(x = all.SNVs$dataset, suffix = 'HGSC')]
}else if(histotype == 'CCOC+EnOC'){
  gr.SNV <- all.SNVs[!endsWith(x = all.SNVs$dataset, suffix = 'HGSC')]
}else{
  gr.SNV <- all.SNVs
}

start(gr.ROI) <- start(gr.ROI) - flanking
end(gr.ROI) <- end(gr.ROI) + flanking
gr.ROI <- subsetByOverlaps(gr.ROI, gr.SNV)
gr.SNV <- subsetByOverlaps(gr.SNV, gr.ROI)

file <- paste0('results/',title,'.',flanking,'bp.filt.txt')
x <- read.table(file = file, header = T, sep = '\t', quote = '')
x <- x[!(x$hgnc_symbol %in% c('FAM83E', 'SPACA4', 'LYPD5', 'C2orf91','HOXD4','HSH2D',paste0('HOXA',c(1:8,10:13)))),]
x <- x[order(x$p.value),]
hgnc_symbols.1 <- as.character(x$hgnc_symbol[x$p.value <= x$p.value[10]])
x <- x[order(x$nSamples, decreasing = T),]
x <- x[x$p.value <= 0.05,]
hgnc_symbols.2 <- as.character(x$hgnc_symbol[x$nSamples >= x$nSamples[10]])
x <- x[x$p.adjust <= 0.25,]
hgnc_symbols.3 <- as.character(x$hgnc_symbol[x$nSamples >= x$nSamples[10]])
hgnc_symbols <- sort(unique(hgnc_symbols.2))

x <- read.table(file = paste0('results/',title,'.',flanking,'bp.filt.txt'), header = T, sep = '\t', quote = '')
x <- x[x$hgnc_symbol %in% hgnc_symbols,]
x <- x[order(x$p.value),]
x <- x[order(x$nSamples, decreasing = T),]
x$hgnc_symbol <- factor(x = x$hgnc_symbol, levels = x$hgnc_symbol)
hgnc_symbols <- x$hgnc_symbol

pb <- txtProgressBar(min = 0, max = length(hgnc_symbols), style = 3)
n <- foreach(i=1:length(hgnc_symbols), .combine = rbind)%dopar%{
  ensembl_gene_id <- bm.genes$ensembl_gene_id[bm.genes$hgnc_symbol == hgnc_symbols[i]]
  .n <- sapply(X = unique(gr.SNV$id), FUN = function(x){
    .gr.ROI <- gr.ROI[gr.ROI$ensembl_gene_id == ensembl_gene_id]
    .gr.SNV <- gr.SNV[gr.SNV$id == x]
    n <- length(subsetByOverlaps(.gr.SNV, .gr.ROI))
    return(n)
  })
  setTxtProgressBar(pb, i)
  return(.n)
}

n <- n[,colSums(x = n)>0]
#sort columns
for(i in nrow(n):1){
  n <- n[,order(n[i,], decreasing = T)]
}
df <- data.frame(nSamples = x$nSamples, 
                 log2FC = log2(x$nSamples/x$exp),
                 p.value = -log10(x$p.value))
df$log2FC[is.infinite(df$log2FC)] <- NA
df$p.value[df$p.value < -log10(0.05)] <- NA
ha <- rowAnnotation(df = df, na_col = 'white',
                    col = list('nSamples' = colorRamp2(breaks = c(min(x$nSamples),max(x$nSamples)), colors = c('#fff5f0', '#67000d')),
                               'log2FC' = colorRamp2(breaks = c(0,max(df$log2FC, na.rm = T)), colors = c('#f7fcf5', '#00441b')),
                               'p.value' = colorRamp2(breaks = c(0,max(df$p.value, na.rm = T)), colors = c('#fcfbfd', '#3f007d'))), 
                    annotation_legend_param = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 6)))

text <- paste0(round(100*rowSums(x = n>0)/length(unique(gr.SNV$id)),0),'%')
ha_text <- rowAnnotation(text = row_anno_text(text), width = max_text_width(text), annotation_name_gp = gpar(fontsize = 6))
rownames(n) <- paste0(hgnc_symbols,'\t',round(100*rowSums(x = n>0)/169,0),'%')
h <- Heatmap(matrix = n, show_column_names = F, col = brewer.pal(n = 7, name = 'Blues'), name = 'nSNVs', row_title = NA, na_col = 'white', 
             column_title = NA, column_title_side = 'bottom', cluster_rows = F, cluster_columns = F, row_names_side = 'left', show_row_names = T, show_heatmap_legend = T, 
             row_names_gp = gpar(fontsize = 6), heatmap_legend_param = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 6)))
w <- 3.125
pdf(file = 'figures/Figure4f.pdf', width = 1.75*w, height = w, useDingbats = F)
print(h+ha)
dev.off()
