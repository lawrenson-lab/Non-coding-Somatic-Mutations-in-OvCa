rm(list = ls())

#bioconductor packages
library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(biomaRt)
library(rtracklayer)

#CRAN packages
library(MASS)
library(ggplot2)
library(poibin)
library(qqman)
library(scales)
library(ggExtra)
library(ggrepel)

source('Rscripts/RegisterCores.R')
source('Rscripts/utils.R')
source('Rscripts/layer_IO.R')
source('Rscripts/layer_business_ncSNVs_PBD_gBMR.R')
source('Rscripts/layer_presentation.R')
source('Rscripts/layer_business.R')
source('Rscripts/ConsensusREs.R')
source('Rscripts/creagUtils.R')

#metadata
samples <- read.table(file = 'data/PCAWG-13_10 Molecular subtypes and Clinical correlates/pcawg_specimen_histology_August2016_v9.txt', header = T, sep = '\t', quote = '')
samples <- samples[samples$histology_abbreviation == 'Ovary-AdenoCA',]
samples.WGS <- samples[samples$specimen_library_strategy =='WGS',]
samples.RNA <- samples[samples$specimen_library_strategy =='RNA-Seq',]
samples <- merge(x = samples.WGS, y = samples.RNA, by = colnames(samples)[c(-4,-6,-7,-28)], suffixes = c('.WGS','.RNA-Seq'))
aliquots.WGS <- read.table(file = 'data/whitelisted_tumour_aliqouts.20160831.txt', header = T, sep = '\t')
aliquots.WGS <- aliquots.WGS[aliquots.WGS$dcc_project_code %in% samples$project_code,]
samples <- merge(x = samples, y = aliquots.WGS, by.x = colnames(samples)[c(2,6,5,1,26,28)], by.y = colnames(aliquots.WGS)[c(-1,-7)])
aliquots.RNA <- read.table(file = 'data/PCAWG/RNA-seq/rnaseq_metadata.tsv', sep = '\t', header = T, quote = "", comment.char = "")
samples <- merge(x = samples, y = aliquots.RNA, by.x = colnames(samples)[c(1,2,4,30,34,3)], by.y = colnames(aliquots.RNA)[c(2,3,4,5,9,37)], suffixes = c('.WGS','.RNA-Seq'))

#WGS
load(file = 'RData/noncodingSNVs.RData')
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/wgEncodeDacMapabilityConsensusExcludable.bed')) == 0]
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/seq.cov1.ONHG19.bed.gz')) == 0]
gr.SNV <- gr.SNV[gr.SNV$id %in% samples$aliquot_id.WGS]

#RNA-Seq
# rna <- pcawg.readRNASeqData()
# rna.sampleIds <- as.vector(as.matrix(rna[1,-1]))
# rna.geneIds <- as.vector(as.matrix(rna[-1,1]))
# colIndexes <- which(rna.sampleIds %in% samples$`aliquot_id.RNA-Seq`)
# .rna <- matrix(data = as.numeric(as.matrix(rna[-1, colIndexes + 1])), nrow = length(rna.geneIds), ncol = length(colIndexes))
# .rna.sampleIds <- rna.sampleIds[colIndexes]
# .rna.donorIds <- sapply(X = .rna.sampleIds, FUN = function(x){
#   .x <- as.character(samples$icgc_donor_id[samples$`aliquot_id.RNA-Seq` == x])
#   return(.x)
# })
# save(list = c('rna.geneIds', '.rna', '.rna.sampleIds', '.rna.donorIds'), file = 'RData/RNA-seq_PCAWG_Ovary-AdenoCA.RData')
load(file = 'RData/RNA-seq_PCAWG_Ovary-AdenoCA.RData')


getPromoters <- function(mart){
  bm.promoters <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'transcript_start', 'strand'), mart = mart)
  bm.promoters$chromosome_name <- paste0('chr', bm.promoters$chromosome_name)
  colnames(bm.promoters)[2:3] <- c('seqnames', 'start')
  bm.promoters$end <- ifelse(test = bm.promoters$strand == 1, yes = bm.promoters$start + 100, no = bm.promoters$start + 1000)
  bm.promoters$start <- ifelse(test = bm.promoters$strand == 1, yes = bm.promoters$start - 1000, no = bm.promoters$start - 100)
  bm.promoters$strand <- ifelse(test = bm.promoters$strand == 1, yes = '+', no = '-')
  gr.promoters <- GRanges(bm.promoters)
  return(gr.promoters)
}

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

minOverlap <- 4
histotype <- 'HGSOC'
creag.init(minOverlap = minOverlap, histotype = histotype)
gr.mapping <- getOCMapping()

ensembl_gene_ids <- unique(subsetByOverlaps(gr.mapping, gr.SNV)$ensembl_gene_id)
is.first <- T
for(i in 1:length(ensembl_gene_ids)){
  .gr <- gr.mapping[gr.mapping$ensembl_gene_id == ensembl_gene_ids[i]]
  mutant.aliquot_id.WGS <-  unique(subsetByOverlaps(gr.SNV, .gr)$id)
  isMutated <- .rna.sampleIds %in% samples$`aliquot_id.RNA-Seq`[samples$aliquot_id.WGS %in% mutant.aliquot_id.WGS]
  .df <- data.frame(ensembl_gene_id = ensembl_gene_ids[i])
  if(sum(isMutated)>1){
    rowIndex <- which(substr(x = rna.geneIds, start = 0, stop = 15) %in% .gr$ensembl_gene_id)
    df <- data.frame(FPKM = .rna[rowIndex,], isMutated = isMutated)
    ks.result <- ks.test(x = df$FPKM[isMutated], y = df$FPKM[!isMutated], alternative = 'two.sided')
    wilcox.result <- wilcox.test(x = df$FPKM[isMutated], y = df$FPKM[!isMutated], alternative = 'two.sided')
    fligner.result <- fligner.test(x = df$FPKM, g = isMutated)
    
    x1 <- median(df$FPKM[isMutated])
    x2 <- median(df$FPKM[!isMutated])
    z1 <- (x1 - mean(df$FPKM))/sd(df$FPKM)
    z2 <- (x2 - mean(df$FPKM))/sd(df$FPKM)
    .df$x1 <- x1
    .df$x2 <- x2
    .df$sd1 <- sd(df$FPKM[isMutated])
    .df$sd2 <- sd(df$FPKM[!isMutated])
    .df$FC <- x1/x2
    .df$deltaGE <- x1-x2
    .df$deltaGEZ <- z1-z2
    .df$geneMeanGE <- mean(df$FPKM)
    .df$geneMedianGE <- median(df$FPKM)
    .df$ks.pvalue <- ks.result$p.value
    .df$wilcox.pvalue <- wilcox.result$p.value
    .df$fligner.pvalue <- fligner.result$p.value
    .df$nMutated <- sum(isMutated)
    write.table(x = .df, file = 'results/pcawgRNAseq_mutantNotMutant_global_new_aim2_v3.txt', append = !is.first, quote = F, sep = '\t', col.names = is.first, row.names = F)
    is.first <- F
  }
}