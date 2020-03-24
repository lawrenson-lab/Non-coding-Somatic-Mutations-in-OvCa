rm(list = ls())

library(rtracklayer)
library(ggplot2)
library(cowplot)

source('Rscripts/RegisterCores.R')
source('Rscripts/layer_IO.R')

gr.ccoc.1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/20160213-ClearCell_241_IP.nodup.tagAlign_x_20160213-ClearCell_241_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.ccoc.2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/20160213-ClearCell_511_IP.nodup.tagAlign_x_20160213-ClearCell_511_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.ccoc.3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/ClearCell-3172_IP.nodup.tagAlign_x_ClearCell-3172_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.ccoc.4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/ClearCell-3547_IP.nodup.tagAlign_x_ClearCell-3547_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.ccoc.5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/ClearCell-3588_IP.nodup.tagAlign_x_ClearCell-3588_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.enoc.1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_291_IP.nodup.tagAlign_x_Endometrioid_291_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.enoc.2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_381_IP.nodup.tagAlign_x_Endometrioid_381_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.enoc.3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_589-Tu_IP.nodup.tagAlign_x_Endometrioid_589-Tu_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.enoc.4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_630_IP.nodup.tagAlign_x_Endometrioid_630_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.enoc.5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_703_IP.nodup.tagAlign_x_Endometrioid_703-Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.hgsoc.1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_229_IP.nodup.tagAlign_x_HGSerous_229_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.hgsoc.2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_270_IP.nodup.tagAlign_x_HGSerous_270_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.hgsoc.3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_429_IP.nodup.tagAlign_x_HGSerous_429_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.hgsoc.4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HG_Serous_550_IP.nodup.tagAlign_x_HG_Serous_550_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.hgsoc.5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_561_IP.nodup.tagAlign_x_HGSerous_561_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.moc.1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_230_IP.nodup.tagAlign_x_Mucinous_230_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.moc.2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_380_IP.nodup.tagAlign_x_Mucinous_380_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.moc.3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_464_IP.nodup.tagAlign_x_Mucinous_464_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.moc.4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_652_IP.nodup.tagAlign_x_Mucinous_652_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
gr.moc.5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous-3724_IP.nodup.tagAlign_x_Mucinous_3724_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')

gr.ccoc.1$ID <- 'CCOC-1'
gr.ccoc.2$ID <- 'CCOC-2'
gr.ccoc.3$ID <- 'CCOC-3'
gr.ccoc.4$ID <- 'CCOC-4'
gr.ccoc.5$ID <- 'CCOC-5'
gr.enoc.1$ID <- 'EnOC-1'
gr.enoc.2$ID <- 'EnOC-2'
gr.enoc.3$ID <- 'EnOC-3'
gr.enoc.4$ID <- 'EnOC-4'
gr.enoc.5$ID <- 'EnOC-5'
gr.hgsoc.1$ID <- 'HGSOC-1'
gr.hgsoc.2$ID <- 'HGSOC-2'
gr.hgsoc.3$ID <- 'HGSOC-3'
gr.hgsoc.4$ID <- 'HGSOC-4'
gr.hgsoc.5$ID <- 'HGSOC-5'
gr.moc.1$ID <- 'MOC-1'
gr.moc.2$ID <- 'MOC-2'
gr.moc.3$ID <- 'MOC-3'
gr.moc.4$ID <- 'MOC-4'
gr.moc.5$ID <- 'MOC-5'

gr <- c(gr.ccoc.1, gr.ccoc.2, gr.ccoc.3, gr.ccoc.4, gr.ccoc.5,
        gr.enoc.1, gr.enoc.2, gr.enoc.3, gr.enoc.4, gr.enoc.5,
        gr.hgsoc.1, gr.hgsoc.2, gr.hgsoc.3, gr.hgsoc.4, gr.hgsoc.5,
        gr.moc.1, gr.moc.2, gr.moc.3, gr.moc.4, gr.moc.5)

grl <- split(x = gr, f = gr$ID)



# for(i in 1:20){
#   combinations <- combn(x = 1:20, m = i)
#   if(ncol(combinations) > choose(n = 20, k = 3)){
#     combinations <- combinations[,sample(x = 1:ncol(combinations), size = 1000, replace = F)]
#   }
#   y <- foreach(j=1:ncol(combinations), .combine = rbind)%dopar%{
#     for(k in 1:i){
#       if(k==1){
#         x <- grl[[combinations[1,j]]]
#       }else{
#         x <- c(x,grl[[combinations[k,j]]])
#       }
#     }
#     x <- reduce(x)
#     nPeaks <- length(x)
#     coverage <- sum(width(reduce(x)))
#     return(data.frame(i,j,nPeaks,coverage))
#   }
#   write.table(x = y, file = 'results/nPeaks_vs_nSamples.txt', append = (i > 1), quote = F, sep = '\t', row.names = F, col.names = (i==1))
# }

z <- read.table(file = 'results/nPeaks_vs_nSamples.txt', header = T, sep = '\t', quote = '')
z$i <- factor(z$i)
z$nPeaks <- z$nPeaks/1000
z$coverage <- z$coverage/1000000

p1 <- ggplot(data = z, mapping = aes(x = i, y = nPeaks)) + 
  geom_boxplot(size = 0.25, outlier.size = 0.25) + xlab('Number of samples') + ylab('Thousands\nof peaks') + theme_classic(base_size = 6)

source('Rscripts/layer_business.R')
hg19.len <- get.hg19.len()
hg19.len$V2 <- hg19.len$V2/1000000
z$coveragePercent <-  100*z$coverage/sum(hg19.len$V2)

p2 <- ggplot(data = z, mapping = aes(x = i, y = coveragePercent)) + geom_boxplot(size = 0.25, outlier.size = 0.25) + xlab('Number of samples') + ylab('Genome\ncoverage (%)') + theme_classic(base_size = 6)

w <- 3.125
pdf(file = 'figures/Figure1b.pdf', width = w, height = w/2, useDingbats = F)
p1 <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p1 <- p1 + theme(plot.margin = margin(0,0,0,0,'in'))
p2 <- p2 + theme(plot.margin = margin(0,0,0,0,'in'))
plot_grid(p1,p2,align = 'hv', axis = 'b', nrow = 2, ncol = 1)
dev.off()

