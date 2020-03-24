rm(list = ls())

#bioconductor packages
library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(biomaRt)
library(rtracklayer)

#CRAN packages
library(MASS)
library(ggExtra)
library(ggplot2)
library(poibin)
library(qqman)
library(scales)

source('Rscripts/RegisterCores.R')
source('Rscripts/utils.R')
source('Rscripts/layer_IO.R')
source('Rscripts/layer_business_ncSNVs_PBD_gBMR.R')
source('Rscripts/layer_presentation.R')
source('Rscripts/layer_business.R')

#metadata
samples <- read.table(file = 'data/PCAWG/PCAWG-13_10 Molecular subtypes and Clinical correlates/pcawg_specimen_histology_August2016_v9.txt', header = T, sep = '\t', quote = '')
samples <- samples[samples$histology_abbreviation == 'Ovary-AdenoCA',]
samples.WGS <- samples[samples$specimen_library_strategy =='WGS',]
samples.RNA <- samples[samples$specimen_library_strategy =='RNA-Seq',]
samples <- merge(x = samples.WGS, y = samples.RNA, by = colnames(samples)[c(-4,-6,-7,-28)], suffixes = c('.WGS','.RNA-Seq'))
aliquots.WGS <- read.table(file = 'data/PCAWG/whitelisted_tumour_aliqouts.20160831.txt', header = T, sep = '\t')
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

bm.genes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                  filters = 'chromosome_name',
                  values = c(1:22,'X'),
                  mart = .ensembl)
overallMedian <- median(.rna[substr(x = rna.geneIds, start = 0, stop = 15) %in% bm.genes$ensembl_gene_id[bm$gene_biotype == 'protein_coding'],])

x <- read.table(file = 'results/pcawgRNAseq_mutantNotMutant_global_new_aim1_mapping.txt', header = T, sep = '\t', quote = '')
x <- merge(x,bm.genes,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
x$log10MedianGE <- log10(x$geneMedianGE)
x$log2FC <- log2(x$FC)
x$binedNDonors <- .bincode(x = x$nMutated, breaks = c(1:12,20), right = F, include.lowest = T)
x$binedNDonors[x$binedNDonors == 12] <- "12+"
x$binedNDonors <- factor(x$binedNDonors, levels = c(1:11,"12+"))
x <- x[order(x$nMutated),]

x$negLog10WilcoxPValue <- -log10(x$wilcox.pvalue)
xMax <- max(abs(x$log2FC[is.finite(x$log2FC)]))
pdf(file = 'figures/Figure3d.pdf', width = w, height = w, useDingbats = F)
p <- ggplot(data = x, mapping = aes(x = log2FC, y = negLog10WilcoxPValue))
p <- p + geom_point() + 
  scale_colour_brewer(palette = 'Paired', name = '#mutated samples') +
  scale_x_continuous(labels = comma, name = expression(-log[2](FC)), limits = c(-xMax,xMax)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'black') + ylab(expression(-log[10](p))) +
  theme(legend.position="top") + theme_classic() +
  guides(col = guide_legend(nrow = 1)) + 
  geom_vline(xintercept = c(-1,1), linetype = 2, col = 'black') + 
  geom_rug(sides = 'b', col=rgb(.5,0,0,alpha=.2))
ggMarginal(p = p, margins = 'x', type = "histogram", xparams = list(breaks=seq(floor(-xMax),xMax,1)))
dev.off()  

