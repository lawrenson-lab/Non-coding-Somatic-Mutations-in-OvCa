rm(list = ls())

library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

source('Rscripts/utils.R')
source('Rscripts/mappingUtils.R')
source('Rscripts/creagUtils.R')

creag.init(minOverlap = 3, histotype = 'ALL')
creags <- getCREAGs(r.cutoff = 0.04, pValue.cutoff = 0.05, rPrime.cutoff = 0.5, dMax = 500000)
length(unique(creags$ensembl_gene_id))
length(unique(creags$enhancer.name))
x.nREs <- sapply(X = unique(creags$ensembl_gene_id), FUN = function(x){
  return(length(unique(creags$enhancer.name[creags$ensembl_gene_id == x])))
})
x.nGenes <- sapply(X = unique(creags$enhancer.name), FUN = function(x){
  return(length(unique(creags$ensembl_gene_id[creags$enhancer.name == x])))
})
(100*sum(!is.na(creags$isGenehancerHit) & creags$isGenehancerHit)/length(creags))
median(x.nREs)
range(x.nREs)
median(x.nGenes)
range(x.nGenes)
mean(x.nREs)
mean(x.nGenes)

df1 <- data.frame(ensembl_gene_id = unique(creags$ensembl_gene_id), n = x.nREs)
df2 <- data.frame(enhancer.name = unique(creags$enhancer.name), n = x.nGenes)

pdf('figures/Figure1i.pdf', width = 2, height = 2, pointsize = 12)
ggplot(data = df1, mapping = aes(x = n)) + geom_histogram(binwidth = 1, fill = '#008B8B', col = 'white', size = 0.1) + 
  xlab('Number of associated REs') +
  ylab('Number of genes') +
  scale_y_log10() +
  theme(text = element_text(family = 'sans', size = 6), 
        panel.grid = element_blank(), panel.background = element_blank(), 
        axis.title = element_text(family = 'sans', size = 6),
        axis.text = element_text(family = 'sans', size = 6))
dev.off()

pdf('figures/Figure1j.pdf', width = 2, height = 2, pointsize = 12)
ggplot(data = df2, mapping = aes(x = n)) + geom_histogram(binwidth = 1, fill = '#008B8B', col = 'white', size = 0.1) + 
  xlab('Number of associated genes') +
  ylab('Number of REs') +
  scale_y_log10() +
  theme(text = element_text(family = 'sans', size = 6), 
        panel.grid = element_blank(), panel.background = element_blank(),
        axis.title = element_text(family = 'sans', size = 6),
        axis.text = element_text(family = 'sans', size = 6))
dev.off()
