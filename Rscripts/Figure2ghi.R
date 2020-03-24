rm(list = ls())

library(ggplot2)

x <- read.table(file = 'results/PPP1R3B KD/CCOC qPCR all cell lines.csv', header = T, sep = ',', quote = '')
x$Cell.line.name <- factor(x = x$Cell.line.name, levels = x$Cell.line.name[order(x$Average, decreasing = T)])
colnames(x)[1] <- 'Histotype'

w <- 180/25.4/3
h <- w/2
# pdf(file = 'figures/CCOC qPCR all cell lines.pdf', width = w, height = h, useDingbats = F)
ggplot(data = x, mapping = aes(x = Cell.line.name, y = Average, ymax = Average + SD, ymin = Average - SD, fill = Histotype)) + 
  geom_bar(stat = 'identity', position = 'dodge') + geom_errorbar() + theme_classic() + 
  scale_fill_manual(values = c('#e41a1c', '#ccebc5', '#4daf4a')) +
  xlab('Cell line') + ylab('PPP1R3B expression') +
  theme(axis.title = element_text(size = 6), axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5),
        legend.justification=c(0,1), 
        legend.position=c(0.85, 0.95),
        legend.background = element_blank(),
        legend.key = element_blank(), legend.key.size = unit(1, 'mm'))
# dev.off()

w <- 2*180/25.4/3/3
h <- w
x$Histotype <- factor(x = x$Histotype, levels = c('CCOC', 'HGSOC', 'FTSEC'))
pdf(file = 'figures/Figure2g.pdf', width = w, height = h, useDingbats = F)
ggplot(data = x, mapping = aes(x = Histotype, y = Average, color = Histotype)) + 
  geom_boxplot(outlier.colour = NA, show.legend = F) + geom_jitter(width = 0.25/2/2, show.legend = F, size = 1) +
  scale_color_manual(values = c('#e41a1c', '#4daf4a', '#ff7f00')) + 
  theme_classic() + theme(axis.title = element_text(size = 6), axis.text = element_text(size = 5)) + ylab('PPP1R3B expression')
dev.off()



x <- read.table(file = 'results/PPP1R3B KD/CCOC GLYCOGEN ASSAY.csv', header = T, sep = ',', quote = '')
x$Experiment <- factor(x = x$Experiment, levels = c('Parental', 'NTC', 'PPP1R3B KD'))
pdf(file = 'figures/Figure2i.pdf', width = w, height = h, useDingbats = F)
ggplot(data = x, mapping = aes(x = Cell.line, y = Average, ymax = Average + SD, ymin = Average - SD, fill = Experiment)) + 
  geom_bar(stat = 'identity', position = 'dodge') + geom_errorbar(stat = 'identity', position = 'dodge') +
  xlab('Cell line') + ylab('Relative glycogen level') + 
  scale_fill_manual(values = c('#1f78b4', '#a6cee3', '#33a02c')) +
  theme_classic() + 
  theme(axis.title = element_text(size = 6), axis.text = element_text(size = 5),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5),
        legend.justification=c(0,1), 
        legend.position=c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(), legend.key.size = unit(1, 'mm'))
dev.off()


x <- read.table(file = 'results/PPP1R3B KD/2019-02-27 091019_QuantStudio 12K Flex_export.csv', header = T, sep = ',', quote = '')
x$Experiment <- factor(x = x$Experiment, levels = c('Parental', 'NTC', 'PPP1R3B KD'))
pdf(file = 'figures/Figure2h.pdf', width = w, height = h, useDingbats = F)
ggplot(data = x, mapping = aes(x = Cell.line, y = Average, ymax = Average + SD, ymin = Average - SD, fill = Experiment)) + 
  geom_bar(stat = 'identity', position = 'dodge') + geom_errorbar(stat = 'identity', position = 'dodge') +
  xlab('Cell line') + ylab('Relative PPP1R3B expression') + 
  scale_fill_manual(values = c('#1f78b4', '#a6cee3', '#33a02c')) +
  theme_classic() + 
  theme(axis.title = element_text(size = 6), axis.text = element_text(size = 5),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5),
        legend.justification=c(0,1), 
        legend.position=c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(), legend.key.size = unit(1, 'mm'))
dev.off()
