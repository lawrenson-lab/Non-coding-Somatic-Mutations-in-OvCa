library(DiffBind)
library(foreach)
library(doMC)
library(VariantAnnotation)
registerDoMC(cores = 4)

getLocation <- function(gr){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  promoters <- promoters(txdb, upstream = 1000, downstream = 100, columns=c("gene_id"))
  location <- rep(NA, length(gr))
  bOverlapsPromoters <- countOverlaps(gr, promoters) > 0
  location[bOverlapsPromoters] <- 'promoter'
  location[!bOverlapsPromoters] <- 'enhancer'
  return(location)
}

.getHistotypeSpecificH3K27acSites <- function(histotype, minOverlap = 3, method = DBA_EDGER, th = 0.05, fold = 3){
  load('RData/diffBind_session03_step01.RData')
  if(histotype=='CCOC'){
    peaks <- 1:5
  }else if(histotype=='EnOC'){
    peaks <- 6:10
  }else if(histotype=='HGSOC'){
    peaks <- 11:15
  }else if(histotype=='MOC'){
    peaks <- 16:20
  }
  dba <- dba.peakset(DBA = dba, peaks = peaks, consensus = DBA_TISSUE, minOverlap = minOverlap)
  gr.consensus <- dba.peakset(DBA = dba, peaks = 21, bRetrieve = T)
  load(paste0('RData/diffBind_session03_step04_', tolower(histotype), '.RData'))
  dba.DB <- dba.report(DBA = dba, method = method, th = th, fold = fold)
  dba.DB <- dba.DB[(dba.DB$Fold>0 & countOverlaps(dba.DB, gr.consensus)>0)|(dba.DB$Fold<0 & countOverlaps(dba.DB, gr.consensus)==0)]
  return(dba.DB)
}

getAllHistotypeSpecificH3K27acSites <- function(minOverlap = 3, method = DBA_EDGER, th = 0.05, fold = 3){
  .gr.ccoc <- .getHistotypeSpecificH3K27acSites(histotype = 'CCOC', minOverlap, method, th, fold)
  .gr.enoc <- .getHistotypeSpecificH3K27acSites(histotype = 'EnOC', minOverlap, method, th, fold)
  .gr.hgsc <- .getHistotypeSpecificH3K27acSites(histotype = 'HGSOC', minOverlap, method, th, fold)
  .gr.muoc <- .getHistotypeSpecificH3K27acSites(histotype = 'MOC', minOverlap, method, th, fold)
  
  gr.ccoc <- .gr.ccoc[countOverlaps(.gr.ccoc, c(granges(.gr.enoc), granges(.gr.hgsc), granges(.gr.muoc)))==0]
  gr.enoc <- .gr.enoc[countOverlaps(.gr.enoc, c(granges(.gr.ccoc), granges(.gr.hgsc), granges(.gr.muoc)))==0]
  gr.hgsc <- .gr.hgsc[countOverlaps(.gr.hgsc, c(granges(.gr.ccoc), granges(.gr.enoc), granges(.gr.muoc)))==0]
  gr.muoc <- .gr.muoc[countOverlaps(.gr.muoc, c(granges(.gr.ccoc), granges(.gr.enoc), granges(.gr.hgsc)))==0]
  
  names(mcols(gr.ccoc))[2:3] <- c('Conc_Hist', 'Conc_Not')
  names(mcols(gr.enoc))[2:3] <- c('Conc_Hist', 'Conc_Not')
  names(mcols(gr.hgsc))[2:3] <- c('Conc_Hist', 'Conc_Not')
  names(mcols(gr.muoc))[2:3] <- c('Conc_Hist', 'Conc_Not')
  
  gr.ccoc$histotype <- 'CCOC'
  gr.enoc$histotype <- 'EnOC'
  gr.hgsc$histotype <- 'HGSOC'
  gr.muoc$histotype <- 'MOC'
  
  gr <- c(gr.ccoc, gr.enoc, gr.hgsc, gr.muoc)
  return(gr)
}

getH3K27acTumorConsensus <- function(){
  load('RData/diffBind_session03_step01.RData')
  gr <- dba.peakset(DBA = dba, peaks = 1:20, minOverlap = 20, bRetrieve = T)
  return(gr)
}

filterGenes <- function(x, cpmCutoff = 1, minSamples = 4){
  cols <- 7:ncol(x)
  x <- x[substr(x$Geneid, 0, 15) %in% bm$ensembl_gene_id[bm$gene_biotype=='protein_coding'],]
  x <- x[substr(rownames(x),1,15) %in% bm$ensembl_gene_id[bm$entrezgene %in% mcols(genes)$gene_id],]
  x <- x[rowSums(cpm(x[,cols]) > cpmCutoff) >= minSamples,]
  return(x)
}

filterSamples <- function(x){
  return(x[,-c(26:30)])
}

.getDEGenes <- function(x, histotype, foldChangeCutoff = 1.2, fdrCutoff = 0.1, positiveFold = T, negativeFold = T){
  countData <- x[,7:ncol(x)]
  condition <- factor(startsWith(colnames(countData), histotype))
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
  dds <- DESeq(dds)
  r <- results(object = dds, contrast = c('condition', TRUE, FALSE))
  if(positiveFold & negativeFold){
    deGenes <- (!is.na(r$log2FoldChange)) & (abs(r$log2FoldChange) >= foldChangeCutoff) & (!is.na(r$padj)) & (r$padj <= fdrCutoff)
  }else if(positiveFold & !negativeFold){
    deGenes <- (!is.na(r$log2FoldChange)) & (r$log2FoldChange >= foldChangeCutoff) & (!is.na(r$padj)) & (r$padj <= fdrCutoff)
  }else if(!positiveFold & negativeFold){
    deGenes <- (!is.na(r$log2FoldChange)) & (r$log2FoldChange <= -foldChangeCutoff) & (!is.na(r$padj)) & (r$padj <= fdrCutoff)
  }
  return(deGenes)
}

getDEGenes <- function(x, positiveFold = T, negativeFold = T){
  deGenesCC <- .getDEGenes(x = x, histotype = 'CCOC', positiveFold = positiveFold, negativeFold = negativeFold)
  deGenesEn <- .getDEGenes(x = x, histotype = 'EnOC', positiveFold = positiveFold, negativeFold = negativeFold)
  deGenesHG <- .getDEGenes(x = x, histotype = 'HGSOC', positiveFold = positiveFold, negativeFold = negativeFold)
  deGenesMu <- .getDEGenes(x = x, histotype = 'MOC', positiveFold = positiveFold, negativeFold = negativeFold)
  selectedGenesCC <- which(deGenesCC & !deGenesEn & !deGenesHG & !deGenesMu)
  selectedGenesEn <- which(!deGenesCC & deGenesEn & !deGenesHG & !deGenesMu)
  selectedGenesHG <- which(!deGenesCC & !deGenesEn & deGenesHG & !deGenesMu)
  selectedGenesMu <- which(!deGenesCC & !deGenesEn & !deGenesHG & deGenesMu)
  index <- c(selectedGenesCC, selectedGenesEn, selectedGenesHG, selectedGenesMu)
  histotype <- c(rep('CCOC',length(selectedGenesCC)),
                 rep('EnOC',length(selectedGenesEn)),
                 rep('HGSOC',length(selectedGenesHG)),
                 rep('MOC',length(selectedGenesMu)))
  df <- data.frame(index, histotype)
  return(df)
}

getNormalizedCounts <- function(x, gene_length, scale = T){
  normalizedCounts <- cpm(x, gene.length = gene_length)
  if(scale){
    normalizedCounts <- t(scale(t(normalizedCounts)))
  }
  return(normalizedCounts)
}

getRawCounts <- function(filterGenes = T){
  x <- readRawCounts()
  if(filterGenes){
    x <- filterGenes(x)  
  }
  x <- filterSamples(x)
}

getChIPSeq_RNASeq_cor <- function(.external_gene_name, windowSize = 100000, method = 'pearson'){
  .bm <- bm[bm$external_gene_name == .external_gene_name,]
  .normalizedCounts <- normalizedCounts[substr(rownames(normalizedCounts),0,15) %in% .bm$ensembl_gene_id,]
  .genes <- genes[genes$gene_id %in% .bm$entrezgene]
  start(.genes) <- start(.genes)-windowSize
  end(.genes) <- end(.genes)+windowSize
  .dba.DB <- subsetByOverlaps(dba.DB, .genes)
  .x <- as.matrix(mcols(.dba.DB)[c(1:7,9:20)])
  .y <- as.vector(as.matrix(.normalizedCounts))
  if(length(.y)>19){
    .y <- as.vector(colMeans(as.matrix(.normalizedCounts)))
  }
  if(nrow(.x)==0 | length(.y)==0){
    return(NA)
  }
  .df <- foreach(j=1:nrow(.x), .combine = rbind) %dopar% {
    .cor <- cor.test(.x[j,],.y, method = method)
    .p.value <- .cor$p.value
    .estimate <- .cor$estimate
    return(c(.p.value, .estimate))
  }
  if(length(.dba.DB)==1){
    .dba.DB$p.value <- .df[1]
    .dba.DB$estimate <- .df[2]
  }else{
    .dba.DB$p.value <- .df[,1]
    .dba.DB$estimate <- .df[,2]
  }
  return(.dba.DB)
}

getTopEnhancer <- function(.external_gene_name, windowSize = 100000){
  .dba.DB <- getChIPSeq_RNASeq_cor(.external_gene_name, windowSize)
  if(is.na(.dba.DB)){
    return(NA)
  }
  .i <- which(.dba.DB$estimate == max(.dba.DB$estimate))
  return(.dba.DB[.i])
}

getRawPeaks <- function(){
  load('RData/diffBind_session03_step01.RData')
  manifest <- dba.show(DBA = dba)
  x <- foreach(i=1:nrow(manifest), .combine = rbind) %dopar% {
    gr <- dba.peakset(DBA = dba, peaks = i, bRetrieve = T)
    Histotype <- rep(manifest$Tissue[i],length(gr))
    sampleId <- rep(manifest$ID[i],length(gr))
    peakWidth <- width(gr)
    score <- gr$V8
    .x <- data.frame(Histotype, sampleId, peakWidth, score)
    return(.x)
  }
  return(x)
}

getChIPSeqQC <- function(){
  x <- readChIPSeqQC()
  x$nReadsCtl <- x$nReadsCtl/1000000
  x$nReadsRep1 <- x$nReadsRep1/1000000
  return(x)
}

getRandomGR <- function(gr, keep.seqnames = TRUE){
  seqnames <- factor(x = paste0('chr',c(1:22,'X')), levels = paste0('chr',c(1:22,'X')))
  hg19.len <- read.table(file = 'data/hg19.len')
  hg19.len <- hg19.len[hg19.len$V1 %in% seqnames,]
  hg19.len$V1 <- factor(x = hg19.len$V1, levels = seqnames)
  hg19.len <- hg19.len[order(hg19.len$V1),]
  
  if(keep.seqnames){
    .seqnames <- factor(as.character(seqnames(gr)), levels = levels(seqnames))
  }else{
    .seqnames <- sample(x = seqnames, size = length(x = gr), replace = T, prob = (hg19.len$V2/10)/sum(hg19.len$V2/10))
  }
  gr.random <- foreach(i=1:length(seqnames), .combine = c) %do% {
    .n <- sum(.seqnames == seqnames[i])
    .max <- hg19.len$V2[i]
    .start <- round(x = runif(n = .n, min = 1, max = .max), digits = 0)
    .width <- width(gr[.seqnames == seqnames[i]])
    .gr <- GRanges(seqnames = seqnames[i], ranges = IRanges(start = .start, width = .width))
    return(.gr)
  }
  return(gr.random)
}

getSuperEnhancersByHistotype <- function(histotype = 'CCOC'){
  gr <- NA
  initialWD <- getwd()
  path <- 'data/superEnhancers'
  setwd(path)
  if(histotype == 'CCOC'){
    gr1 <- homer.readSuperEnhancers(filename ='20160213-ClearCell_241_IP.nodup.superEnhancers.txt')
    gr2 <- homer.readSuperEnhancers(filename = '20160213-ClearCell_511_IP.nodup.superEnhancers.txt')
    gr3 <- homer.readSuperEnhancers(filename = 'ClearCell-3172_IP.nodup.superEnhancers.txt')
    gr4 <- homer.readSuperEnhancers(filename = 'ClearCell-3547_IP.nodup.superEnhancers.txt')
    gr5 <- homer.readSuperEnhancers(filename = 'ClearCell-3588_IP.nodup.superEnhancers.txt')
  }else if(histotype == 'EnOC'){
    gr1 <- homer.readSuperEnhancers(filename ='Endometrioid_291_IP.nodup.superEnhancers.txt')
    gr2 <- homer.readSuperEnhancers(filename = 'Endometrioid_381_IP.nodup.superEnhancers.txt')
    gr3 <- homer.readSuperEnhancers(filename = 'Endometrioid_589-Tu_IP.nodup.superEnhancers.txt')
    gr4 <- homer.readSuperEnhancers(filename = 'Endometrioid_630_IP.nodup.superEnhancers.txt')
    gr5 <- homer.readSuperEnhancers(filename = 'Endometrioid_703_IP.nodup.superEnhancers.txt')
  }else if(histotype == 'HGSOC'){
    gr1 <- homer.readSuperEnhancers(filename ='HGSerous_229_IP.nodup.superEnhancers.txt')
    gr2 <- homer.readSuperEnhancers(filename = 'HGSerous_270_IP.nodup.superEnhancers.txt')
    gr3 <- homer.readSuperEnhancers(filename = 'HGSerous_429_IP.nodup.superEnhancers.txt')
    gr4 <- homer.readSuperEnhancers(filename = 'HG_Serous_550_IP.nodup.superEnhancers.txt')
    gr5 <- homer.readSuperEnhancers(filename = 'HGSerous_561_IP.nodup.superEnhancers.txt')
  }else if(histotype == 'MOC'){
    gr1 <- homer.readSuperEnhancers(filename ='Mucinous_230_IP.nodup.superEnhancers.txt')
    gr2 <- homer.readSuperEnhancers(filename = 'Mucinous_380_IP.nodup.superEnhancers.txt')
    gr3 <- homer.readSuperEnhancers(filename = 'Mucinous_464_IP.nodup.superEnhancers.txt')
    gr4 <- homer.readSuperEnhancers(filename = 'Mucinous_652_IP.nodup.superEnhancers.txt')
    gr5 <- homer.readSuperEnhancers(filename = 'Mucinous-3724_IP.nodup.superEnhancers.txt')
  }
  gr <- c(gr1, gr2, gr3, gr4, gr5)
  setwd(initialWD)
  return(gr)
}

getLocalMutationRate <- function(gr.H3K27ac, gr.snv, windowSize = 1000000){
  .gr.H3K27ac <- gr.H3K27ac
  start1 <- start(.gr.H3K27ac) - windowSize
  end1 <- start(.gr.H3K27ac) - 1
  start2 <- end(.gr.H3K27ac) + 1
  end2 <- end(.gr.H3K27ac) + windowSize
  start(.gr.H3K27ac) <- start1
  end(.gr.H3K27ac) <- end1
  width1 <- width(.gr.H3K27ac)
  .x1 <- countOverlaps(.gr.H3K27ac, gr.snv)
  end(.gr.H3K27ac) <- end2
  start(.gr.H3K27ac) <- start2
  width2 <- width(.gr.H3K27ac)
  .x2 <- countOverlaps(.gr.H3K27ac, gr.snv)
  .x <- 1000*(.x1+.x2)/(width1+width2)
  return(.x)
}

getLocalMutationRate2 <- function(gr.H3K27ac, gr.snv, windowSize = 1000000){
  .gr.H3K27ac <- gr.H3K27ac
  start1 <- start(.gr.H3K27ac) - windowSize
  end1 <- start(.gr.H3K27ac) - 1
  start2 <- end(.gr.H3K27ac) + 1
  end2 <- end(.gr.H3K27ac) + windowSize
  start(.gr.H3K27ac) <- start1
  end(.gr.H3K27ac) <- end1
  width1 <- width(.gr.H3K27ac)
  .x1 <- countOverlaps(gr.snv, .gr.H3K27ac)$id
  end(.gr.H3K27ac) <- end2
  start(.gr.H3K27ac) <- start2
  width2 <- width(.gr.H3K27ac)
  .x2 <- countOverlaps(gr.snv, .gr.H3K27ac)$id
  .x <- 1000*(length(unique(c(.x1,.x2))))/(width1+width2)
  return(.x)
}

getDips <- function(.gr.narrowPeaks, gr.bigWig, order = 250){
  .name <- .gr.narrowPeaks$name
  .gr.bigWig <- subsetByOverlaps(query = gr.bigWig, subject = .gr.narrowPeaks)
  z <- foreach(i=1:length(.gr.bigWig), .combine = rbind)%do%{
    .s <- start(.gr.bigWig[i]) 
    .w <- width(.gr.bigWig[i])
    .x <- seq(.s, .s+.w-1)
    .y <- rep(.gr.bigWig$score[i],.w)
    .z <- cbind(.x,.y)
    return(.z)
  }
  
  y1 <- z[,2]
  y2 <- ma(x = y1, order = order, centre = F)
  y2 <- ma(x = y2, order = order, centre = F)
  y3 <- c(NA,y2[-1] - y2[-length(y2)])
  iLocalMin <- which(y3[-1] > 0 & y3[-length(y3)] < 0)
  iLocalMax <- which(y3[-1] < 0 & y3[-length(y3)] > 0)
  
  gr.dips <- foreach(i=1:length(iLocalMin), .combine = c, .errorhandling = "remove")%do%{
    iMin <- iLocalMin[i]
    iMax1 <- max(iLocalMax[(iLocalMin[i] - iLocalMax) > 0])
    iMax2 <- min(iLocalMax[(iLocalMax - iLocalMin[i]) > 0])
    xMin <- z[iMin,1]
    xMax1 <- z[iMax1,1]
    xMax2 <- z[iMax2,1]
    yMin <- z[iMin,2]
    yMax1 <- z[iMax1,2]
    yMax2 <- z[iMax2,2]
    w <- iMax2 - iMax1
    w1 <- iMin - iMax1
    w2 <- iMax2 - iMin
    d1 <- yMax1 - yMin
    d2 <- yMax2 - yMin
    a1 <- w1*d1/2
    a2 <- w2*d2/2
    a3 <- w*min(d1,d2) + w*(max(d1,d2)-min(d1,d2))/2 - a1 - a2
    .gr <- GRanges(seqnames = seqnames(.gr.narrowPeaks), ranges = IRanges(start = xMax1, end = xMax2), 
                   xMin = xMin, yMax1 = yMax1, yMax2 = yMax2, yMin = yMin, area = a3, 
                   peakMax = max(z[,2]), peakName = .name)
    return(.gr)
  }
  return(gr.dips)
}

mapEnahcner2Gene <- function(seqnames, start, end, windowSize = 1000000, method = 'BH', mu = log(3)/10000){
  # mu: decay rate
  gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end))
  chipSeq <- mcols(subsetByOverlaps(query = dba.DB, subject = gr))
  chipSeq <- chipSeq[,-8]
  
  .gr <- gr
  start(.gr) <- start(.gr)-windowSize
  end(.gr) <- end(.gr)+windowSize
  
  .genes <- bm.tss
  .genes <- .genes[.genes$ensembl_gene_id %in% substr(x = rownames(normalizedCounts), start = 0, stop = 15)]
  .genes <- subsetByOverlaps(query = .genes, subject = .gr)
  if(length(.genes) > 0){
    c <- foreach(i=1:length(.genes), .combine = c)%do%{
      .gene <- .genes[i]
      if(overlapsAny(query = gr, subject = .gene)){
        .gene$d <- 0
      }else{
        .gene$d <- min(abs(start(gr) - start(.gene)),abs(end(gr) - start(.gene)))
      }
      row <- substr(x = rownames(normalizedCounts), start = 0, stop = 15) %in% .gene$ensembl_gene_id
      .cor <- cor.test(x = as.matrix(chipSeq)[1,], y = normalizedCounts[row,])
      .gene$r <- .cor$estimate
      .gene$p.value <- .cor$p.value
      return(.gene)
    }
    d <- c$d
    w <- 2*exp(-mu*d)/(1+exp(-mu*d))
    rPrime <- w*c$r
    p.adjust <- p.adjust(p = c$p.value, method = method)
    
    c$w <- w
    c$rPrime <- rPrime
    c$p.adjust <- p.adjust
    c$enhancer <- rep(paste0(seqnames,':',start,'-',end),length(c))
    c$gene <- paste0(seqnames(c),':',start(c),'-',end(c))
    df <- as.data.frame(c)
  }else{
    df <- NA
  }
  return(df)
}

get.hg19.len <- function(){
  seqnames <- factor(x = paste0('chr',c(1:22,'X')), levels = paste0('chr',c(1:22,'X')))
  hg19.len <- read.table(file = 'data/hg19.len')
  hg19.len <- hg19.len[hg19.len$V1 %in% seqnames,]
  hg19.len$V1 <- factor(x = hg19.len$V1, levels = seqnames)
  hg19.len <- hg19.len[order(hg19.len$V1),]
  return(hg19.len)
}

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