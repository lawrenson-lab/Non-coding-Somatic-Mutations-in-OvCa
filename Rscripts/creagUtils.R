library(biomaRt)
library(ggplot2)
library(rtracklayer)

source('Rscripts/RegisterCores.R')

creag.init <- function(minOverlap, histotype){
  gr.tads <<- getTADs(extend = TRUE)
  .ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = 'grch37.ensembl.org')
  .ensembl <- useDataset('hsapiens_gene_ensembl',mart=.ensembl)
  bm.genes <<- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'chromosome_name',
                     values = c(1:22,'X'),
                     mart = .ensembl)
  bm.tss <<- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'transcription_start_site'),
                   filters = 'chromosome_name',
                   values = c(1:22,'X'),
                   mart = .ensembl)
  
  if(histotype %in% c('CCOC', 'EnOC', 'HGSOC', 'MOC')){
    gr.consensusHGSOC <<- getConsensusREs(minOverlap = minOverlap, histotype = histotype)
  }else if(histotype == 'CCOC+EnOC'){
    gr.consensusHGSOC <<- c(getConsensusREs(minOverlap = minOverlap, histotype = 'CCOC'),
                            getConsensusREs(minOverlap = minOverlap, histotype = 'EnOC'))
    gr.consensusHGSOC <<- reduce(gr.consensusHGSOC)
  }else if(histotype == 'ALL'){
    gr.consensusHGSOC <<- c(getConsensusREs(minOverlap = minOverlap, histotype = 'CCOC'),
    getConsensusREs(minOverlap = minOverlap, histotype = 'EnOC'),
    getConsensusREs(minOverlap = minOverlap, histotype = 'HGSOC'),
    getConsensusREs(minOverlap = minOverlap, histotype = 'MOC'))
    gr.consensusHGSOC <<- reduce(gr.consensusHGSOC)
  }
  start(gr.consensusHGSOC) <<- start(gr.consensusHGSOC)
  end(gr.consensusHGSOC) <<- end(gr.consensusHGSOC)
  load(paste0('results/enhancerGeneMapDB_spearman_genehancerHits.Rdata'))
  df.mapping$enhancer.name <- paste0(df.mapping$seqnames, ':', df.mapping$start, '-', df.mapping$end)
  df.mapping <<- df.mapping
}

plotCREAG <- function(hgnc_symbol, cor.method, rCutoff = 0.5, pValueCutoff = 0.05, rPrimeCutoff = 0){
  df <- read.table(file = paste0('results/enhnacerGeneMappingDB.txt'), header = F, sep = '\t', quote = '')
  colnames(df) <- c('seqnames', 'start', 'end', 'ensembl_gene_id', 'mu', 'cor.method', 'd', 'r', 'p.value', 'w', 'rPrime')
  df <- df[df$cor.method == cor.method,]
  df$enhancer.name <- paste0(df$seqnames, ':', df$start, '-', df$end)
  
  gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$start, end = df$end))
  gr <- subsetByOverlaps(query = gr, subject = gr.consensusHGSOC)
  df <- df[df$enhancer.name %in% paste0(seqnames(gr), ':', start(gr), '-', end(gr)),]
  
  source(file = 'Rscripts/manuallyCurated_enhancerNames.R')
  df <- setEnhancerName(df)
  
  .df <- df[df$ensembl_gene_id == bm.genes$ensembl_gene_id[bm.genes$hgnc_symbol == hgnc_symbol] & 
              df$r >= rCutoff & 
              df$p.value <= pValueCutoff &
              df$rPrime >= rPrimeCutoff,]
  .df <- .df[order(.df$r, decreasing = T),]
  df <- df[df$enhancer.name %in% .df$enhancer.name,]
  df$enhancer.name <- factor(x = df$enhancer.name, levels = unique(df$enhancer.name))
  df$ensembl_gene_id <- factor(x = df$ensembl_gene_id, levels = unique(bm.tss$ensembl_gene_id[order(bm.tss$transcription_start_site)]))
  df <- df[order(df$ensembl_gene_id),]
  df$r[df$p.value > pValueCutoff | df$r < 0] <- NA
  df$r <- factor(x = round(df$r, 1), levels = seq(0,1,0.1))
  df <- merge(x = df, y = bm.genes, all.x = T, sort = F)
  df$hgnc_symbol <- factor(x = df$hgnc_symbol, levels = unique(df$hgnc_symbol))
  p <- ggplot(df, aes(hgnc_symbol, enhancer.name)) + 
    geom_tile(aes(fill = r), colour = "white") + 
    scale_fill_brewer(palette = 'Blues') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    xlab('Gene\nsorted by genomic location') + ylab('Regulatory elements\nsorted by genomic location') +
    ggtitle(label = paste0(hgnc_symbol, ' CREAG'), subtitle = paste0('cor.method = ', cor.method, ', r > ',rCutoff,', p < ',pValueCutoff, ', rPrime > ',rPrimeCutoff))
  return(p)
}

getCREAG <- function(hgnc_symbol, r.cutoff = 0.4, pValue.cutoff = 0.05, FDR = 1, rPrime.cutoff = 0.5, dMax = 500000){
  creag <- df.mapping
  mu <- log(3)/dMax
  w <- 2*exp(-mu*creag$d)/(1+exp(-mu*creag$d))
  creag$rPrime <- w*creag$r
  
  creag <- creag[creag$ensembl_gene_id == bm.genes$ensembl_gene_id[bm.genes$hgnc_symbol == hgnc_symbol],]
  creag$p.adjust <- p.adjust(p = creag$p.value, method = 'fdr')
  creag <- creag[creag$r >= r.cutoff &
                   creag$p.value <= pValue.cutoff & 
                   creag$p.adjust <= FDR &
                   creag$rPrime >= rPrime.cutoff,]
  gr.mapping <- GRanges(creag)
  gr.mapping <- subsetByOverlaps(query = gr.mapping, subject = gr.consensusHGSOC)
  return(gr.mapping)
}

getCREAGs <- function(r.cutoff = 0.4, pValue.cutoff = 0.05, rPrime.cutoff = 0.5, dMax = 500000){
  creag <- df.mapping
  mu <- log(3)/dMax
  w <- 2*exp(-mu*creag$d)/(1+exp(-mu*creag$d))
  creag$rPrime <- w*creag$r
  creag <- creag[creag$r >= r.cutoff &
                   creag$p.value <= pValue.cutoff & 
                   creag$rPrime >= rPrime.cutoff,]
  gr.mapping <- GRanges(creag)
  genome(gr.mapping) <- 'hg19'
  genome(gr.consensusHGSOC) <- 'hg19'
  gr.mapping <- subsetByOverlaps(gr.mapping, gr.consensusHGSOC)
  return(gr.mapping)
}
