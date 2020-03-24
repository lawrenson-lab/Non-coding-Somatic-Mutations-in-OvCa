rm(list = ls())
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(poibin)

source('Rscripts/layer_IO.R')
source('Rscripts/layer_business.R')

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cds <- cds(txdb)

cellLine <- 'MCF7'

gr.TFBS <- import.bed(con = paste0('data/remap2018_',cellLine,'_all_macs2_hg19_v1_2.bed.gz'))
gr.TFBS <- gr.TFBS[seqnames(gr.TFBS) %in% paste0('chr',c(1:22,'X'))]
TF.ChIPseq <- unique(gr.TFBS$name)
histologyAbbreviations <- as.character(pcawg.getHistologyAbbreviations())
histologyAbbreviations <- histologyAbbreviations[histologyAbbreviations != 'Bone-Leiomyo']

minOverlap <- 5
histotype <- 'HGSOC'
flanking <- 500
gr.ROI <- getConsensusREs(minOverlap = minOverlap, histotype = histotype)
start(gr.ROI) <- start(gr.ROI) - flanking
end(gr.ROI) <- end(gr.ROI) + flanking
nGlobal <- sum(width(reduce(gr.ROI)))

# for(col in 1:length(histologyAbbreviations)){
#   print(paste0(col, ': ', histologyAbbreviations[col]))
#   gr.SNV <- pcawg.getSNVsByHistologyAbbreviation(histology_abbreviation = histologyAbbreviations[col])
#   gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = cds) == 0]
#   gr.SNV <- subsetByOverlaps(gr.SNV, gr.ROI)
#   save(list = 'gr.SNV', file = paste0('RData/noncodingSNVs_PCAWG_',histologyAbbreviations[col],'_ROI.RData'))
# }


is.first <- T
for(col in 1:length(histologyAbbreviations)){
  print(paste0(col, ': ', histologyAbbreviations[col]))
  load(file = paste0('RData/noncodingSNVs_PCAWG_',histologyAbbreviations[col],'_ROI.RData'))
  
  p <- sapply(X = unique(gr.SNV$aliquot_id), FUN = function(x){
    .gr.SNV <- gr.SNV[gr.SNV$aliquot_id == x]
    pi <- length(.gr.SNV)/nGlobal
    return(pi)
  })
  for(row in 1:length(TF.ChIPseq)){
    x <- gr.TFBS[gr.TFBS$name == TF.ChIPseq[row]]
    ni <- sum(width(reduce(x)))
    x <- subsetByOverlaps(x, gr.ROI, type = 'within')
    nj <- sum(width(reduce(x)))
    .gr.SNV <- subsetByOverlaps(gr.SNV, x)
    nSNVs <- length(.gr.SNV)
    nSamples <- length(unique(.gr.SNV$aliquot_id))
    pj <- p
    pp <- 1-(1-pj)^nj
    .k <-  nSamples
    pvalue.poibin <- 1-ppoibin(kk = .k-1, pp = pp, method = 'DFT-CF')
    pvalue.binom <- 1-pbinom(q = .k-1, size = length(pp), prob = mean(pp))
    pvalue.pois <- 1-ppois(q = .k-1, lambda = length(pp)*mean(pp))
    expVal.poibin <- qpoibin(qq = c(0.5), pp = pp)
    expVal.binom <- qbinom(prob = mean(pp), size = length(pp), p = c(0.5))
    expVal.pois <- qpois(p = c(0.5), lambda = length(pp)*mean(pp))
    sigVal.poibin <- qpoibin(qq = c(0.95), pp = pp)
    sigVal.binom <- qbinom(prob = mean(pp), size = length(pp), p = c(0.95))
    sigVal.pois <- qpois(p = c(0.95), lambda = length(pp)*mean(pp))
    df <- data.frame(cancerType = histologyAbbreviations[col], 
                     TF.ChIPseq = TF.ChIPseq[row], 
                     totalTFBS = ni, 
                     totalSamples = length(pp), 
                     bmr = mean(pp),
                     activeTFBS = nj, 
                     nSNVs, nSamples, 
                     pvalue.poibin, pvalue.binom, pvalue.pois, 
                     expVal.poibin, expVal.binom, expVal.pois,
                     sigVal.poibin, sigVal.binom, sigVal.pois)
    write.table(x = df, file = paste0('results/enrichmentTFBS_',cellLine,'.txt'), append = !is.first, quote = F, sep = '\t', col.names = is.first, row.names = F)
    is.first <- F
  }
}