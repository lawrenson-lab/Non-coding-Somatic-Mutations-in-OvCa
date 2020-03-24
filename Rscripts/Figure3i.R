rm(list = ls())
library(pzfx)
library(ggplot2)
tables <- pzfx_tables(path = 'chr6Enhacer/Chr6_KO_Single & Bulk Cell.pzfx')
x <- read_pzfx(path = 'chr6Enhacer/Chr6_KO_Single & Bulk Cell.pzfx', table = 'Chr6_SHIN3_Single Cell')
gene <- x[,1]
nGenes <- length(gene)
x <- as.matrix(x[,-1])

df <- data.frame(gene,
                 clone = c(rep('OR1C1 KO', 3*nGenes), rep('Enh WT', 4*nGenes), rep('Enh PKO', 6*nGenes), rep('Enh CKO', 3*nGenes)),
                 biologicalReplicate = c(rep(1:3,each = nGenes), rep(1:4, each = nGenes), rep(1:6, each = nGenes), rep(1:3, each = nGenes)),
                 rbind(x[,1:2], x[,3:4], x[,5:6], x[,13:14], x[,15:16], x[,17:18], x[,19:20], x[,25:26], x[,27:28], x[,29:30], x[,31:32], x[,33:34], x[,35:36], x[,37:38], x[,39:40], x[,41:42]))

colnames(df)[4:5] <- c('relativeExpression.1', 'relativeExpression.2')
df$relativeExpression <- rowMeans(df[,4:5])
df$sd <- apply(df[,4:6], 1, sd, na.rm = TRUE)
df$gene <- factor(x = df$gene, levels = c('ZSCAN16', 'ZSCAN12', 'HIST1H2AI', 'ZKSCAN3', 'ZSCAN31'))
df$clone <- factor(x = df$clone, levels = c('OR1C1 KO', 'Enh WT', 'Enh PKO', 'Enh CKO'))
w <- 180/25.4/2
pdf(file = 'figures/Figure3i.pdf', width = w, height = w/2, useDingbats = F)
ggplot(data = df, mapping = aes(x = gene, y = relativeExpression, fill = clone)) + 
  geom_boxplot(outlier.colour = NA, show.legend = F) +
  geom_jitter(shape = 21, position = position_dodge(width = 0.75)) + 
  geom_hline(yintercept = 1, lty = 2) +
  scale_fill_brewer(palette = 'Paired') +
  ylab('Relative Expression') +
  theme_classic() + theme(axis.title.x = element_blank(), legend.title = element_blank(), 
                          text = element_text(size = 6), 
                          line = element_line(colour = 'black', size = 1), 
                          rect = element_rect(colour = 'black', size = 1))
dev.off()
