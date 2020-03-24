rm(list = ls())

#bioconductor packages
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(rtracklayer)

#CRAN packages
library(poibin)
library(qqman)

source('Rscripts/RegisterCores.R')
source('Rscripts/layer_IO.R')
source('Rscripts/layer_business_ncSNVs_PBD_gBMR.R')
source('Rscripts/layer_business.R')
source(file = 'Rscripts/ConsensusREs.R')

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cds <- cds(txdb)

load(file = 'RData/noncodingSNVs.RData')
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/wgEncodeDacMapabilityConsensusExcludable.bed')) == 0]
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/seq.cov1.ONHG19.bed.gz')) == 0]
all.SNVs <- gr.SNV

# parameters
histotypes <- c('CCOC', 'EnOC', 'HGSOC', 'CCOC+EnOC','CCOC+EnOC+HGSOC')
minOverlaps <- 1:5
flankings <- seq(from=0, to = 500, by = 50)

for(histotype in histotypes){
  for(minOverlap in minOverlaps){
    for(flanking in flankings){
      print(paste(histotype, minOverlap, flanking))
      if(histotype %in% c('CCOC','EnOC','HGSOC')){
        gr.ROI <- getConsensusREs(minOverlap = minOverlap, histotype = histotype)
        if(histotype == 'CCOC'){
          gr.SNV <- all.SNVs[endsWith(x = all.SNVs$dataset, suffix = 'CCOC')]
        }else if(histotype == 'EnOC'){
          gr.SNV <- all.SNVs[endsWith(x = all.SNVs$dataset, suffix = 'ENOC')]
        }else{
          gr.SNV <- all.SNVs[endsWith(x = all.SNVs$dataset, suffix = 'HGSC')]
        }
      }else if(histotype == 'CCOC+EnOC'){
        gr.ROI <- reduce(c(getConsensusREs(minOverlap = minOverlap, histotype = 'CCOC'),
                           getConsensusREs(minOverlap = minOverlap, histotype = 'EnOC')))
        gr.SNV <- all.SNVs[!endsWith(x = all.SNVs$dataset, suffix = 'HGSC')]
      }else{
        gr.ROI <- reduce(c(getConsensusREs(minOverlap = minOverlap, histotype = 'CCOC'),
                           getConsensusREs(minOverlap = minOverlap, histotype = 'EnOC'),
                           getConsensusREs(minOverlap = minOverlap, histotype = 'HGSOC')))
        gr.SNV <- all.SNVs
      }
      start(gr.ROI) <- start(gr.ROI) - flanking
      end(gr.ROI) <- end(gr.ROI) + flanking
      
      n <- sum(width(reduce(gr.ROI)))
      gr.ROI <- subsetByOverlaps(gr.ROI, gr.SNV)
      gr.SNV <- subsetByOverlaps(gr.SNV, gr.ROI)
      
      #get background mutation rate
      p <- sapply(X = unique(gr.SNV$id), FUN = function(x){
        .gr.SNV <- gr.SNV[gr.SNV$id == x]
        pi <- length(subsetByOverlaps(.gr.SNV, gr.ROI))/n
        return(pi)
      })
      mcols <- foreach(i=1:length(gr.ROI), .combine = rbind)%dopar%{
        x <- gr.ROI[i]
        nj <- sum(width(reduce(x)))
        .gr.SNV <- subsetByOverlaps(gr.SNV, x)
        nSNVs <- length(.gr.SNV)
        nSamples <- length(unique(.gr.SNV$id))
        pj <- p
        pp <- 1-(1-pj)^nj
        .k <-  nSamples
        .pvalue <- 1-ppoibin(kk = .k-1, pp = pp, method = 'DFT-CF')
        return(c(nj,nSNVs,nSamples,.pvalue))
      }
      colnames(mcols) <- c('nj', 'nSNV', 'nSamples', 'p.value')
      mcols(gr.ROI) <- mcols
      gr.ROI$p.adjust <- p.adjust(p = gr.ROI$p.value, method = 'fdr')
      write.table(x = gr.ROI, file = paste0('results/aim1/aim1_FMREs_',histotype,'_minOverlap',minOverlap,'_flanklin',flanking,'bp.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)
    }
  }
}