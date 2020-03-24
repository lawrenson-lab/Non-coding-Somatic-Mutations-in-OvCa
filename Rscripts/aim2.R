rm(list = ls())

library(poibin)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
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

histotypes <- c('HGSOC', 'CCOC', 'EnOC', 'CCOC+EnOC', 'ALL')
minOverlap <- 4
flanking <- 0
for(histotype in histotypes){
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
  
  n <- sum(width(reduce(gr.ROI, ignore.strand = T)))
  gr.ROI <- subsetByOverlaps(gr.ROI, gr.SNV)
  gr.SNV <- subsetByOverlaps(gr.SNV, gr.ROI)
  
  p <- sapply(X = unique(gr.SNV$id), FUN = function(x){
    .gr.SNV <- gr.SNV[gr.SNV$id == x]
    pi <- length(subsetByOverlaps(.gr.SNV, gr.ROI))/n
    return(pi)
  })
  
  mcols <- foreach(i=1:length(ensembl_gene_ids), .combine = rbind)%dopar%{
    x <- gr.ROI[gr.ROI$ensembl_gene_id == ensembl_gene_ids[i]]
    nj <- sum(width(reduce(x, ignore.strand = T)))
    .gr.SNV <- subsetByOverlaps(gr.SNV, x)
    nSNVs <- length(.gr.SNV)
    nSamples <- length(unique(.gr.SNV$id))
    pj <- p
    pp <- 1-(1-pj)^nj
    .k <-  nSamples
    .pvalue <- 1-ppoibin(kk = .k-1, pp = pp, method = 'DFT-CF')
    exp <- qpoibin(qq = 0.5, pp = pp)
    result <- c(as.character(ensembl_gene_ids[i]), nj,nSNVs,nSamples,.pvalue, exp, length(pp), length(reduce(x, ignore.strand = T)), length(reduce(x[strand(x)=='*'])), length(reduce(x[strand(x)!='*'])))
    return(result)
  }
  mcols <- as.data.frame(mcols)
  colnames(mcols) <- c('ensembl_gene_id', 'nj', 'nSNV', 'nSamples', 'p.value', 'exp', 'nTotSamples', 'nREs', 'nNonPromoters', 'nPromoters')
  
  bm.genes <- getBM(attributes=c('ensembl_gene_id', 'description', 'gene_biotype', 'hgnc_symbol', 
                                 'chromosome_name', 'start_position', 'end_position', 'strand'),
                    filters = 'chromosome_name',
                    values = c(1:22,'X'),
                    mart = .ensembl)
  
  mcols <- merge(x = bm.genes, y = mcols, by = 'ensembl_gene_id', all.y = T)
  mcols$p.value <- as.numeric(as.character(mcols$p.value))
  mcols$nSamples <- as.numeric(as.character(mcols$nSamples))
  mcols$nNonPromoters <- as.numeric(as.character(mcols$nNonPromoters))
  mcols2 <- mcols[mcols$nSamples > 0 & mcols$nNonPromoters > 0,]
  mcols2$p.adjust <- p.adjust(p = mcols2$p.value, method = 'fdr')
  write.table(x = mcols2, file = paste0('results/',title,'.',flanking,'bp.filt.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}
