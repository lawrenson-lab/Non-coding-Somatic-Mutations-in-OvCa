library(rtracklayer)

source('Rscripts/layer_IO.R')

getConsensusREs <- function(minOverlap = 3, histotype){
  if(histotype == 'CCOC'){
    gr1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/20160213-ClearCell_241_IP.nodup.tagAlign_x_20160213-ClearCell_241_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/20160213-ClearCell_511_IP.nodup.tagAlign_x_20160213-ClearCell_511_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/ClearCell-3172_IP.nodup.tagAlign_x_ClearCell-3172_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/ClearCell-3547_IP.nodup.tagAlign_x_ClearCell-3547_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/ClearCell-3588_IP.nodup.tagAlign_x_ClearCell-3588_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
  }else if(histotype == 'EnOC'){
    gr1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_291_IP.nodup.tagAlign_x_Endometrioid_291_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_381_IP.nodup.tagAlign_x_Endometrioid_381_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_589-Tu_IP.nodup.tagAlign_x_Endometrioid_589-Tu_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_630_IP.nodup.tagAlign_x_Endometrioid_630_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Endometrioid_703_IP.nodup.tagAlign_x_Endometrioid_703-Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
  }else if(histotype == 'HGSOC'){
    gr1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_229_IP.nodup.tagAlign_x_HGSerous_229_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_270_IP.nodup.tagAlign_x_HGSerous_270_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_429_IP.nodup.tagAlign_x_HGSerous_429_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HG_Serous_550_IP.nodup.tagAlign_x_HG_Serous_550_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/HGSerous_561_IP.nodup.tagAlign_x_HGSerous_561_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
  }else if(histotype == 'MOC'){
    gr1 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_230_IP.nodup.tagAlign_x_Mucinous_230_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr2 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_380_IP.nodup.tagAlign_x_Mucinous_380_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr3 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_464_IP.nodup.tagAlign_x_Mucinous_464_INPUT.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr4 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous_652_IP.nodup.tagAlign_x_Mucinous_652_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
    gr5 <- read.narrowPeak(filename = 'chipseq_hg19/out/peak/macs2/overlap/Mucinous-3724_IP.nodup.tagAlign_x_Mucinous_3724_Input.nodup.tagAlign.naive_overlap.filt.narrowPeak.gz')
  }else{
    return(NULL)
  }
  gr.union <- c(gr1, gr2, gr3, gr4, gr5)
  gr.union <- gr.union[seqnames(gr.union) %in% paste0('chr',c(1:22,'X'))]
  seqlevels(gr.union) <- paste0('chr',c(1:22,'X'))
  grl <- split(x = gr.union, f = seqnames(gr.union))
  gr.consensus <- unlist(GRangesList(lapply(X = grl, FUN = function(x){
    dividers <- unique(sort(c(start(x),end(x))))
    start <- dividers[-length(dividers)]
    end <- dividers[-1]
    x <- GRanges(seqnames = unique(seqnames(x)), ranges = IRanges(start = start, end = end), strand = '*')
    x$overlaps <- rowSums(cbind(
      countOverlaps(query = x, subject = gr1, type = 'within') > 0,
      countOverlaps(query = x, subject = gr2, type = 'within') > 0,
      countOverlaps(query = x, subject = gr3, type = 'within') > 0,
      countOverlaps(query = x, subject = gr4, type = 'within') > 0,
      countOverlaps(query = x, subject = gr5, type = 'within') > 0))
    return(reduce(x[x$overlaps >= minOverlap]))
  })))
  return(gr.consensus)
}