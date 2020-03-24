rm(list = ls())

library(ggplot2)
library(ggrepel)
library(scales)


path <- 'results/RoseRes'
dirs <- list.dirs(path = path, recursive = F)
selectedPlots <- c(grep(pattern = 'ClearCell_241', x = dirs),
                   grep(pattern = 'ClearCell_511', x = dirs),
                   grep(pattern = 'ClearCell_3172', x = dirs),
                   grep(pattern = 'ClearCell_3547', x = dirs),
                   grep(pattern = 'ClearCell_3588', x = dirs),
                   grep(pattern = 'Endometrioid_291', x = dirs),
                   grep(pattern = 'Endometrioid_381', x = dirs),
                   grep(pattern = 'Endometrioid_589_Tu', x = dirs),
                   grep(pattern = 'Endometrioid_630', x = dirs),
                   grep(pattern = 'Endometrioid_703', x = dirs),
                   grep(pattern = 'HGSerous_229', x = dirs),
                   grep(pattern = 'HGSerous_270', x = dirs),
                   grep(pattern = 'HGSerous_429', x = dirs),
                   grep(pattern = 'HGSerous_550', x = dirs),
                   grep(pattern = 'HGSerous_561', x = dirs),
                   grep(pattern = 'Mucinous_230', x = dirs),
                   grep(pattern = 'Mucinous_380', x = dirs),
                   grep(pattern = 'Mucinous_464', x = dirs),
                   grep(pattern = 'Mucinous_652', x = dirs),
                   grep(pattern = 'Mucinous_3724', x = dirs))
plot.title <- c(paste0('CCOC-',1:5),
                paste0('EnOC-',1:5),
                paste0('HGSOC-',1:5),
                paste0('MOC-',1:5))
plot.list <- list()
plotNumber <- 1
for(i in selectedPlots){
  title <- strsplit(x = dirs[i], split = '/')[[1]][3]
  print(title)
  file <- paste(dirs[i],list.files(path = dirs[i], pattern = "*_AllEnhancers.table.txt"), sep = '/')
  x <- read.table(file = file, header = T)
  file <- paste(dirs[i],list.files(path = dirs[i], pattern = "*_SuperEnhancers_ENHANCER_TO_GENE.ivyAnnotation.txt"), sep = '/')
  x.SEs <- read.table(file = file, header = T, sep = '\t')
  
  x$label <- NA
  selectedGenes <- c('MUC16', 'CRABP2', 'EPCAM', 'ESR1', 'SPON1', 'WFDC2', 'HNF1A', 'HNF1B', 'IGF2', 'CDH6', 'MKI67', 'KISS1', 'STI4', 'MSLN', 'MIF', 'MMP7', 'CDKN1A', 'TP53', 'PAX8', 'PGR', 'SLPI', 'TACSTD2', 'WT1')
  selectedSEs <- sapply(X = selectedGenes, FUN = function(.x){
    selectedSE <- c(which(as.character(x.SEs$hgnc_symbol) == .x))
    return(min(selectedSE))
  })
  selectedSEs <- sort(unique(selectedSEs[!is.infinite(selectedSEs)]))
  x$label[selectedSEs] <- paste0(selectedSEs, '. ', as.character(x.SEs$hgnc_symbol[selectedSEs]))
  
  x$y <- x[,7] - x[,8]
  x$isSuper <- factor(x$isSuper)
  x$enhancerRank <- length(x$enhancerRank) - x$enhancerRank
  nSuperEnhancers <- sum(x$isSuper == 1)
  yCutoff <- max(x$y[x$isSuper == 0])
  plot.list[[plotNumber]] <- ggplot(data = x, mapping = aes(x = enhancerRank, y = y, col = isSuper, label = label)) + 
    ggtitle(label = plot.title[plotNumber]) +
    geom_line(show.legend = F) + 
    xlab(label = 'Enhancer rank') + 
    ylab(label = 'Enhancer signal (total rpm)') + 
    geom_hline(yintercept = yCutoff, linetype = 2, col = 'grey', size = 0.2) + 
    geom_vline(xintercept = nrow(x) - nSuperEnhancers, linetype = 2, col = 'grey', size = 0.2) + 
    scale_color_brewer(palette = 'Paired') + 
    scale_y_continuous(labels = comma) +
    geom_text_repel(show.legend = F, 
                    col = 'black', 
                    direction    = "y",
                    segment.size = 0.2, 
                    nudge_x = -nrow(x)/2, 
                    hjust = 0.5,
                    size = 2) +
    theme(text = element_text(size = 7),
          axis.title = element_text(size = 7), 
          axis.text.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.background = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = NA), 
          panel.grid = element_blank(),
          plot.title = element_text(size = 7))
  plotNumber <- plotNumber+1
}

library(cowplot)
pdf(file = 'figures/FigureS4a.pdf', width = 10, height = 7.5, compress = T)
plot_grid(plotlist = plot.list, align = 'hv', axis = 'btrl', nrow = 4, ncol = 5)
dev.off()



