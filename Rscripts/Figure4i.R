rm(list = ls())

library(poibin)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(VariantAnnotation)
library(biomaRt)
library(qqman)

source('Rscripts/RegisterCores.R')
source('Rscripts/creagUtils.R')
source('Rscripts/utils.R')
source('Rscripts/ConsensusREs.R')
source('Rscripts/layer_presentation.R')

getPValue <- function(nj, nSamples, p){
  pj <- p
  pp <- 1-(1-pj)^nj
  .k <-  nSamples
  .pvalue <- 1-ppoibin(kk = .k-1, pp = pp, method = 'DFT-CF')
  return(.pvalue)
}

getExpected <- function(nj, nSamples, p){
  pj <- p
  pp <- 1-(1-pj)^nj
  .k <-  nSamples
  expectedValue <- qpoibin(qq = 0.5, pp = pp)
  return(expectedValue)
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

#read active regions
minOverlap <- 4
histotype <- 'HGSOC'
flanking <- 0
gr.ROI <- getConsensusREs(minOverlap = minOverlap, histotype = histotype)
start(gr.ROI) <- start(gr.ROI) - flanking
end(gr.ROI) <- end(gr.ROI) + flanking
nGlobal <- sum(width(reduce(gr.ROI)))

histology_abbreviations <- pcawg.getHistologyAbbreviations()
histology_abbreviations <- histology_abbreviations[histology_abbreviations != 'Bone-Leiomyo']
i <- which(histology_abbreviations == 'Ovary-AdenoCA')
histology_abbreviation <- as.character(histology_abbreviations[i])
aliquots <- pcawg.getAliquotsByHistologyAbbreviation(histology_abbreviation = histology_abbreviation)
aliquot_ids <- sort(unique(aliquots$aliquot_id))

#read SNVs
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cds <- cds(txdb)
load('RData/noncodingSNVs_PCAWG_Ovary-AdenoCA.RData')
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/wgEncodeDacMapabilityConsensusExcludable.bed')) == 0]
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = import.bed(con = 'data/seq.cov1.ONHG19.bed.gz')) == 0]
gr.SNV <- gr.SNV[countOverlaps(query = gr.SNV, subject = cds) == 0]


#read CNV
load(file = 'RData/CNV_PCAWG_Ovary-AdenoCA.RData')
df.cnv.tmp <- df.cnv[-1,]
colnames(df.cnv.tmp) <- as.character(t(df.cnv[1,]))
df.cnv <- df.cnv.tmp
rm(list = 'df.cnv.tmp')
df.cnv <- df.cnv[df.cnv[,1] == 'PAX8',]

#calculate background probability
p <- sapply(X = aliquot_ids, FUN = function(x){
  .gr.SNV <- gr.SNV[gr.SNV$aliquot_id == x]
  pi <- length(subsetByOverlaps(.gr.SNV, gr.ROI))/nGlobal
  return(pi)
})

#count TEAD4 BS mutations per sample
cellLine <- 'MCF7'
gr.TFBS <- import.bed(con = paste0('data/remap2018_',cellLine,'_all_macs2_hg19_v1_2.bed.gz'))
gr.TFBS <- gr.TFBS[seqnames(gr.TFBS) %in% paste0('chr',c(1:22,'X'))]
TF.ChIPseq <- unique(gr.TFBS$name)
row <- grep(pattern = 'TEAD4', x = TF.ChIPseq)
gr.TFBS <- gr.TFBS[gr.TFBS$name == TF.ChIPseq[row]]
.gr.ROI <- subsetByOverlaps(gr.TFBS, gr.ROI, type = 'within')
n.TEAD4.BS <- sapply(X = aliquot_ids, FUN = function(x){
  .gr.SNV <- gr.SNV[gr.SNV$aliquot_id == x]
  .n <- length(subsetByOverlaps(.gr.SNV, .gr.ROI))
  return(.n)
})
observed.TEAD4.BS <- sum(n.TEAD4.BS>0)
pvalue.TEAD4.BS <- getPValue(nj = sum(width(reduce(.gr.ROI))), nSamples = observed.TEAD4.BS, p = p)
expected.TEAD4.BS <- getExpected(nj = sum(width(reduce(.gr.ROI))), nSamples = observed.TEAD4.BS, p = p)
.gr.ROI.TEAD4 <- .gr.ROI

#count PAX8 BS mutations per sample
gr.TFBS <- import.bed(con = 'data/PAX8_ChIP-seq.bed')
gr.TFBS <- gr.TFBS[gr.TFBS$name == 'Kuromachi']
.gr.ROI <- subsetByOverlaps(gr.TFBS, gr.ROI, type = 'within')
n.PAX8.BS <- sapply(X = aliquot_ids, FUN = function(x){
  .gr.SNV <- gr.SNV[gr.SNV$aliquot_id == x]
  .n <- length(subsetByOverlaps(.gr.SNV, .gr.ROI))
  return(.n)
})
observed.PAX8.BS <- sum(n.PAX8.BS>0)
pvalue.PAX8.BS <- getPValue(nj = sum(width(reduce(.gr.ROI))), nSamples = observed.PAX8.BS, p = p)
expected.PAX8.BS <- getExpected(nj = sum(width(reduce(.gr.ROI))), nSamples = observed.PAX8.BS, p = p)
.gr.ROI.PAX8 <- .gr.ROI

creag.init(minOverlap = minOverlap, histotype = histotype)
.gr.ROI <- getOCMapping()
hgnc_symbol <- 'PAX8'
ensembl_gene_id <- bm.genes$ensembl_gene_id[bm.genes$hgnc_symbol == hgnc_symbol]
.gr.ROI <- .gr.ROI[.gr.ROI$ensembl_gene_id == ensembl_gene_id]
n.PAX8.CREAG <- sapply(X = aliquot_ids, FUN = function(x){
  .gr.SNV <- gr.SNV[gr.SNV$aliquot_id == x]
  .n <- length(subsetByOverlaps(.gr.SNV, .gr.ROI))
  return(.n)
})
observed.PAX8.CREAG <- sum(n.PAX8.CREAG>0)
pvalue.PAX8.CREAG <- getPValue(nj = sum(width(reduce(.gr.ROI, ignore.strand = T))), nSamples = observed.PAX8.CREAG, p = p)
expected.PAX8.CREAG <- getExpected(nj = sum(width(reduce(.gr.ROI, ignore.strand = T))), nSamples = observed.PAX8.CREAG, p = p)

data1 <- rbind(as.numeric(t(df.cnv[,4:ncol(df.cnv)])), n.PAX8.CREAG, n.TEAD4.BS, n.PAX8.BS)
n <- matrix(data = as.numeric(data1), nrow = 4, byrow = F)
cnv <- as.numeric(t(df.cnv[,4:ncol(df.cnv)]))
rownames(n) <- paste0(c('CNV', 'PAX8 BS','TEAD4 BS', 'PAX8 CREAG'),'\t',round(100*rowSums(x = n!=0)/ncol(n),0),'%')
colnames(n) <- aliquot_ids
#sort columns
for(i in nrow(n):1){
  order <- order(n[i,], decreasing = T)
  n <- n[,order]
  cnv <- cnv[order]
}
col_ha <- HeatmapAnnotation(df = data.frame(cnv = n[1,]), col = list(cnv = c("2" =  "#ca0020",
                                                                             "1" =  "#f4a582",
                                                                             "0" =  "#f7f7f7",
                                                                             "-1" =  "#92c5de",
                                                                             "-2" =  "#0571b0")),
                            show_legend = F)
h <- Heatmap(matrix = n[-1,], show_column_names = F, col = colorRamp2(breaks = c(0,2), colors = c('#f7f7f7', '#0571b0')), name = 'nSNVs', row_title = NA, na_col = 'white',
             column_title = NA, column_title_side = 'bottom', cluster_rows = F, cluster_columns = F, row_names_side = 'left', show_row_names = T, show_heatmap_legend = F,
             row_names_gp = gpar(fontsize = 6), heatmap_legend_param = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 6)), top_annotation = col_ha)

nSamples <- rowSums(n[-1,]>0)
df <- data.frame(nSamples = nSamples)

pdf(file = 'figures/Figure4i.pdf', width = 4, height = 1)
print(h)
dev.off()
