rm(list = ls())
library(reshape2)
library(ggplot2)
cellLine <- 'MCF7'
cancerType <- 'Ovary-AdenoCA'
df <- read.table(file = paste0('results/enrichmentTFBS_',cellLine,'.txt'), header = T, sep = '\t', quote = '')
df$expID <- matrix(data = unlist(strsplit(x = as.character(df$TF.ChIPseq), split = '.', fixed = T)), ncol = 3, byrow = T)[,1]
df$TF <- matrix(data = unlist(strsplit(x = as.character(df$TF.ChIPseq), split = '.', fixed = T)), ncol = 3, byrow = T)[,2]
df$cellType <- matrix(data = unlist(strsplit(x = as.character(df$TF.ChIPseq), split = '.', fixed = T)), ncol = 3, byrow = T)[,3]
df <- df[df$cellType == cellLine,]
df <- df[df$cancerType == cancerType,]
for(TF in unique(df$TF)){
  if(sum(df$TF == TF)>1){
    i <- which(df$TF == TF & df$pvalue.poibin > min(df$pvalue.poibin[df$TF == TF]))
    df <- df[-i,]
  }
}
df$FDR <- p.adjust(p = df$pvalue.poibin, method = 'fdr')
fdr <- 0.05
.x1 <- -log10(min(df$pvalue.poibin[df$FDR > fdr]))/2
.x2 <- -log10(max(df$pvalue.poibin[df$FDR <= fdr]))/2
df$pvalue.poibin <- -log10(df$pvalue.poibin)


library(ggrepel)

df$FC <- df$nSamples/df$expVal.poibin
df$TF <- factor(x = df$TF, levels = df$TF[order(df$FC, decreasing = T)])
df <- df[!is.nan(df$FC),]
pdf(file = 'figures/Figure3j.pdf', width = 188/25.4/3, height = 188/25.4/3)
ggplot(data = df, mapping = aes(x = FC, y = pvalue.poibin, label = TF)) + 
  geom_vline(xintercept = 1.3, linetype = 3) +
  geom_point() + 
  geom_text_repel(data = df[df$FC > 1.3 & df$pvalue.poibin >= -log10(0.05),], size = 2) + 
  geom_hline(yintercept = -log10(0.05), linetype = 3) + 
  xlab('Fold change (obs/exp)') +
  ylab(expression(-log[10](p))) + 
  theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 6), axis.text = element_text(size = 6))
dev.off()

