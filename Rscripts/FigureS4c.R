rm(list = ls())

library(edgeR)
library(matrixStats)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(reshape2)
library(ggplot2)

source('Rscripts/RegisterCores.R')
source('Rscripts/utils.R')
source('Rscripts/layer_IO.R')
source('Rscripts/layer_business.R')

filterGenes <- function(x, cpmCutoff = 1, minSamples = 4){
  cols <- 7:ncol(x)
  x <- x[rowSums(cpm(x[,cols]) > cpmCutoff) >= minSamples,]
  return(x)
}

ensembl.hg38 = useMart("ensembl",dataset="hsapiens_gene_ensembl")
bm.genes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                  filters = 'chromosome_name',
                  values = c(1:22,'X'),
                  mart = ensembl.hg38)

bm.tss <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand', 'transcript_biotype'),
                filters = 'chromosome_name',
                values = c(1:22,'X'),
                mart = .ensembl)
transcript_biotype <- sort(unique(bm.tss$transcript_biotype))
transcript_biotype <- transcript_biotype[-grep(pattern = 'pseudogene', x = transcript_biotype)]
bm.tss$chromosome_name <- paste0('chr', bm.tss$chromosome_name)
bm.tss$strand <- ifelse(test = bm.tss$strand == 1, yes = '+', no = '-')
colnames(bm.tss) <- c('ensembl_gene_id', 'seqnames', 'start', 'end', 'strand', 'transcript_biotype')
gr.tss <- GRanges(bm.tss)
gr.tss <- gr.tss[gr.tss$transcript_biotype %in% transcript_biotype]

x <- readRawCounts()
x <- filterSamples(x)
x <- filterGenes(x)
normalizedCounts <- getNormalizedCounts(x[,7:ncol(x)], x$Length, F)
medianCPM <- rowMedians(x = normalizedCounts, na.rm = T)
df.medianCPM <- data.frame(ensembl_gene_id = substr(rownames(normalizedCounts), 0, 15), medianCPM)
df.medianCPM <- df.medianCPM[order(medianCPM, decreasing = T),]
df.medianCPM <- merge(x = df.medianCPM, y = bm.genes, sort = F, all.x = T)


path <- 'results/metascapeSEs/'
files <- list.files(path = path, pattern = '*.txt', all.files = T)
for(i in 1:length(files)){
  title <- substr(x = files[i], start = 0, stop = nchar(files[i])-4)
  print(title)
  file <- paste0(path, files[i])
  x <- read.table(file = file, header = F, as.is = T, col.names = 'x')
  .normalizedCounts <- normalizedCounts[substr(rownames(normalizedCounts), 0, 15) %in% x$x,]
  .normalizedCounts <- t(scale(t(.normalizedCounts)))
  ha <- HeatmapAnnotation(df = data.frame(histology = substr(x = colnames(.normalizedCounts), start = 0, stop = nchar(colnames(.normalizedCounts))-2)), which = 'column', 
                          col = list(histology = c('CCOC' = '#e41a1c', 'EnOC' = '#377eb8', 'HGSOC' = '#4daf4a', 'MOC' = '#984ea3')))
  h <- Heatmap(matrix = .normalizedCounts, clustering_distance_rows = 'euclidean', cluster_columns = F,
               col = colorRamp2(c(-4, 0, 4), c("green", "white", "red")), show_row_names = F, name = paste0(title,'\nGene expression\n(z-score)'), show_row_dend = F, bottom_annotation = ha)
  print(h)
}

selectedFiles <- c(1,2,4,5)
x.union <- sort(unique(foreach(i=selectedFiles, .combine = c)%do%{
  title <- substr(x = files[i], start = 0, stop = nchar(files[i])-4)
  print(title)
  file <- paste0(path, files[i])
  x <- read.table(file = file, header = F, as.is = T, col.names = 'x')
  return(x$x)
}))

x.ccoc <- read.table(file = paste0(path, files[1]), header = F, as.is = T, col.names = 'x')$x
x.enoc <- read.table(file = paste0(path, files[2]), header = F, as.is = T, col.names = 'x')$x
x.hgsoc <- read.table(file = paste0(path, files[4]), header = F, as.is = T, col.names = 'x')$x
x.moc <- read.table(file = paste0(path, files[5]), header = F, as.is = T, col.names = 'x')$x
x.union <- c(x.ccoc, x.enoc, x.hgsoc, x.moc)



.normalizedCounts <- normalizedCounts[match(x = x.union, table = substr(rownames(normalizedCounts), 0, 15)),]
.normalizedCounts <- t(scale(t(.normalizedCounts)))

data <- data.frame(SEAG_histology = c(rep('CCOC', length(x.ccoc)), rep('EnOC', length(x.enoc)), rep('HGSOC', length(x.hgsoc)), rep('MOC', length(x.moc))), .normalizedCounts)
data <- melt(data)
data$sampleHistology <- substr(x = data$variable, start = 0, stop = nchar(as.character(data$variable))-2)

pdf(file = 'figures/FigureS4c.pdf', width = 188/25.4/2, height = 188/25.4/6, useDingbats = F)
ggplot(data = data, mapping = aes(x = SEAG_histology, y = value, col = sampleHistology)) + geom_boxplot(outlier.size = 0.5) + scale_color_brewer(palette = 'Set1') +
  ylab('Gene expression (z-score)') + xlab('SEAG specificity') +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic(base_size = 6) +
  theme(axis.text = element_text(colour = 'black'), legend.position = 'none', plot.margin = margin(0,0,0,0,'in'))
dev.off()

