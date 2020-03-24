getROI <- function(minOverlapRE = 3, histotypes = c('CCOC','EnOC','HGSOC')){
  # ROI
  load('RData/diffBind_session03_step01.RData')
  gr.H3K27ac.ccoc <- dba.peakset(DBA = dba, peaks = 1:5, bRetrieve = T, minOverlap = minOverlapRE)
  gr.H3K27ac.enoc <- dba.peakset(DBA = dba, peaks = 6:10, bRetrieve = T, minOverlap = minOverlapRE)
  gr.H3K27ac.hgsc <- dba.peakset(DBA = dba, peaks = 11:15, bRetrieve = T, minOverlap = minOverlapRE)
  gr.H3K27ac.muoc <- dba.peakset(DBA = dba, peaks = 16:20, bRetrieve = T, minOverlap = minOverlapRE)
  gr.H3K27ac.ccoc <- gr.H3K27ac.ccoc[seqnames(gr.H3K27ac.ccoc) %in% paste0('chr',c(1:22,'X'))]
  gr.H3K27ac.enoc <- gr.H3K27ac.enoc[seqnames(gr.H3K27ac.enoc) %in% paste0('chr',c(1:22,'X'))]
  gr.H3K27ac.hgsc <- gr.H3K27ac.hgsc[seqnames(gr.H3K27ac.hgsc) %in% paste0('chr',c(1:22,'X'))]
  gr.H3K27ac.muoc <- gr.H3K27ac.muoc[seqnames(gr.H3K27ac.muoc) %in% paste0('chr',c(1:22,'X'))]
  load('RData/diffBind_session03_step02.RData')
  dba <- dba.count(DBA = dba, peaks = NULL, score = DBA_SCORE_RPKM_FOLD)
  gr.H3K27ac.all <- dba.peakset(DBA = dba, peaks = 1:20, bRetrieve = T, minOverlap = 1)
  .gr <- NA
  if('CCOC' %in% histotypes){
    if(is.na(.gr)){
      .gr <- granges(gr.H3K27ac.ccoc)
    }else{
      .gr <- c(.gr, granges(gr.H3K27ac.ccoc))
    }
  }
  if('EnOC' %in% histotypes){
    if(is.na(.gr)){
      .gr <- granges(gr.H3K27ac.enoc)
    }else{
      .gr <- c(.gr, granges(gr.H3K27ac.enoc))
    }
  }
  if('HGSOC' %in% histotypes){
    if(is.na(.gr)){
      .gr <- granges(gr.H3K27ac.hgsc)
    }else{
      .gr <- c(.gr, granges(gr.H3K27ac.hgsc))
    }
  }
  if('MOC' %in% histotypes){
    if(is.na(.gr)){
      .gr <- granges(gr.H3K27ac.muoc)
    }else{
      .gr <- c(.gr, granges(gr.H3K27ac.muoc))
    }
  }
  gr.H3K27ac.all <- subsetByOverlaps(query = gr.H3K27ac.all, subject = .gr)
  gr.H3K27ac.all <- gr.H3K27ac.all[seqnames(gr.H3K27ac.all) %in% paste0('chr',c(1:22,'X'))]
  gr.H3K27ac.all$elementType <- getElementType(gr.H3K27ac.all, minOverlapSE = 3)
  return(gr.H3K27ac.all)
}

getROI.PAX8 <- function(reduce = F){
  gr.PAX8.HeyA8 <- read.IDROutput(filename = 'data/HeyA8_PAX8.conservative.IDR0.05.filt.narrowPeak')
  gr.PAX8.IGROV1 <- read.IDROutput(filename = 'data/IGROV_PAX8.conservative.IDR0.05.filt.narrowPeak')
  gr.PAX8.FT33 <- import.bed.withCommentedLines(file = 'data/GSM2107936_FT33_PAX8_peaks.bed')
  gr.PAX8.FT194 <- import.bed.withCommentedLines(file = 'data/GSM2107938_FT194_PAX8_peaks.bed')
  gr.PAX8.FT246 <- import.bed.withCommentedLines(file = 'data/GSM2107940_FT246_PAX8_peaks.bed')
  gr.PAX8.Kuromachi <- import.bed.withCommentedLines(file = 'data/GSM2107942_Kuromachi_PAX8_peaks.bed')
  gr.PAX8.OVSAHO <- import.bed.withCommentedLines(file = 'data/GSM2107945_OVSAHO_PAX8_peaks.bed')
  gr.PAX8.JHOS4 <- import.bed.withCommentedLines(file = 'data/GSM2107947_JHOS4_PAX8_peaks.bed')
  gr.PAX8 <- c(granges(gr.PAX8.HeyA8), granges(gr.PAX8.IGROV1), granges(gr.PAX8.FT33), 
               granges(gr.PAX8.FT194), granges(gr.PAX8.FT246), granges(gr.PAX8.Kuromachi), 
               granges(gr.PAX8.OVSAHO), granges(gr.PAX8.JHOS4))
  if(reduce){
    gr.PAX8 <- reduce(gr.PAX8)
  }
  gr.PAX8$elementType <- getElementType(gr.PAX8, minOverlapSE = 3)
  return(gr.PAX8)
}

getElementType <- function(gr, minOverlapSE = 3){
  gr.se <- getEOCSuperEnhancers(minOverlapSE = minOverlapSE)
  gr.cpgIslandExt <- readCpGIslands()
  tss <- promoters(txdb, upstream = 1000, downstream = 100, columns=c("gene_id"))
  promoters$overlapCGI <- countOverlaps(query = tss, subject = gr.cpgIslandExt) > 0
  
  overlapSE <- countOverlaps(query = gr, subject = gr.se)>0
  overlapPr <- countOverlaps(query = gr, subject = promoters)>0
  overlapCGI <- countOverlaps(query = gr, subject = promoters[promoters$overlapCGI])>0
  
  elementType <- rep(NA, length(gr))
  elementType[overlapSE] <- 'overlapSE'
  elementType[!overlapSE & !overlapPr] <- 'nonPromoters'
  elementType[!overlapSE & overlapPr & overlapCGI] <- 'promotersCGI'
  elementType[!overlapSE & overlapPr & !overlapCGI] <- 'promotersNonCGI'
  elementType <- factor(elementType)
  return(elementType)
}

getEOCSuperEnhancers <- function(minOverlapSE){
  # super enhancers
  gr.se.ccoc.1 <- homer.readSuperEnhancers('data/superEnhancers/20160213-ClearCell_241_IP.nodup.superEnhancers.txt')
  gr.se.ccoc.2 <- homer.readSuperEnhancers('data/superEnhancers/20160213-ClearCell_511_IP.nodup.superEnhancers.txt')
  gr.se.ccoc.3 <- homer.readSuperEnhancers('data/superEnhancers/ClearCell-3172_IP.nodup.superEnhancers.txt')
  gr.se.ccoc.4 <- homer.readSuperEnhancers('data/superEnhancers/ClearCell-3547_IP.nodup.superEnhancers.txt')
  gr.se.ccoc.5 <- homer.readSuperEnhancers('data/superEnhancers/ClearCell-3588_IP.nodup.superEnhancers.txt')
  gr.se.enoc.1 <- homer.readSuperEnhancers('data/superEnhancers/Endometrioid_291_IP.nodup.superEnhancers.txt')
  gr.se.enoc.2 <- homer.readSuperEnhancers('data/superEnhancers/Endometrioid_381_IP.nodup.superEnhancers.txt')
  gr.se.enoc.3 <- homer.readSuperEnhancers('data/superEnhancers/Endometrioid_589-Tu_IP.nodup.superEnhancers.txt')
  gr.se.enoc.4 <- homer.readSuperEnhancers('data/superEnhancers/Endometrioid_630_IP.nodup.superEnhancers.txt')
  gr.se.enoc.5 <- homer.readSuperEnhancers('data/superEnhancers/Endometrioid_703_IP.nodup.superEnhancers.txt')
  gr.se.hgsc.1 <- homer.readSuperEnhancers('data/superEnhancers/HGSerous_229_IP.nodup.superEnhancers.txt')
  gr.se.hgsc.2 <- homer.readSuperEnhancers('data/superEnhancers/HGSerous_270_IP.nodup.superEnhancers.txt')
  gr.se.hgsc.3 <- homer.readSuperEnhancers('data/superEnhancers/HGSerous_429_IP.nodup.superEnhancers.txt')
  gr.se.hgsc.4 <- homer.readSuperEnhancers('data/superEnhancers/HG_Serous_550_IP.nodup.superEnhancers.txt')
  gr.se.hgsc.5 <- homer.readSuperEnhancers('data/superEnhancers/HGSerous_561_IP.nodup.superEnhancers.txt')
  
  gr.se.ccoc <- reduce(c(gr.se.ccoc.1, gr.se.ccoc.2, gr.se.ccoc.3, gr.se.ccoc.4, gr.se.ccoc.5))
  .overlaps.ccoc <- cbind(countOverlaps(query = gr.se.ccoc, subject = gr.se.ccoc.1)>0,
                          countOverlaps(query = gr.se.ccoc, subject = gr.se.ccoc.2)>0,
                          countOverlaps(query = gr.se.ccoc, subject = gr.se.ccoc.3)>0,
                          countOverlaps(query = gr.se.ccoc, subject = gr.se.ccoc.4)>0,
                          countOverlaps(query = gr.se.ccoc, subject = gr.se.ccoc.5)>0)
  
  gr.se.enoc <- reduce(c(gr.se.enoc.1, gr.se.enoc.2, gr.se.enoc.3, gr.se.enoc.4, gr.se.enoc.5))
  .overlaps.enoc <- cbind(countOverlaps(query = gr.se.enoc, subject = gr.se.enoc.1)>0,
                          countOverlaps(query = gr.se.enoc, subject = gr.se.enoc.2)>0,
                          countOverlaps(query = gr.se.enoc, subject = gr.se.enoc.3)>0,
                          countOverlaps(query = gr.se.enoc, subject = gr.se.enoc.4)>0,
                          countOverlaps(query = gr.se.enoc, subject = gr.se.enoc.5)>0)
  
  gr.se.hgsc <- reduce(c(gr.se.hgsc.1, gr.se.hgsc.2, gr.se.hgsc.3, gr.se.hgsc.4, gr.se.hgsc.5))
  .overlaps.hgsc <- cbind(countOverlaps(query = gr.se.hgsc, subject = gr.se.hgsc.1)>0,
                          countOverlaps(query = gr.se.hgsc, subject = gr.se.hgsc.2)>0,
                          countOverlaps(query = gr.se.hgsc, subject = gr.se.hgsc.3)>0,
                          countOverlaps(query = gr.se.hgsc, subject = gr.se.hgsc.4)>0,
                          countOverlaps(query = gr.se.hgsc, subject = gr.se.hgsc.5)>0)
  
  gr.se.ccoc <- gr.se.ccoc[rowSums(.overlaps.ccoc) >= minOverlapSE]
  gr.se.enoc <- gr.se.enoc[rowSums(.overlaps.enoc) >= minOverlapSE]
  gr.se.hgsc <- gr.se.hgsc[rowSums(.overlaps.hgsc) >= minOverlapSE]
  gr.se <- reduce(c(gr.se.ccoc, gr.se.enoc, gr.se.hgsc))
  return(gr.se)
}

getSNVs <- function(filterCoding = TRUE, filterSpliceSites = F){
  gr.snv.pcawg <- pcawg.getSNVsByHistologyAbbreviation(histology_abbreviation = 'Ovary-AdenoCA')
  gr.snv.ccoc <- shah.getSNVsByProject(project = 'CCOC')
  gr.snv.enoc <- shah.getSNVsByProject(project = 'ENOC', remove.EnOC.outlier = T)
  gr.snv.hgsc <- shah.getSNVsByProject(project = 'HGSC')
  gr.snv.all <- c(granges(gr.snv.pcawg), granges(gr.snv.ccoc), granges(gr.snv.enoc), granges(gr.snv.hgsc))
  gr.snv.all$id <- c(as.character(gr.snv.pcawg$aliquot_id),
                     as.character(gr.snv.ccoc$case_id),
                     as.character(gr.snv.enoc$case_id),
                     as.character(gr.snv.hgsc$case_id))
  gr.snv.all$dataset <- c(rep('P.HGSC',length(gr.snv.pcawg)),
                          rep('S.CCOC',length(gr.snv.ccoc)),
                          rep('S.ENOC',length(gr.snv.enoc)),
                          rep('S.HGSC',length(gr.snv.hgsc)))
  gr.snv.all$REF <- c(as.character(gr.snv.pcawg$REF), 
                      as.character(gr.snv.ccoc$ref), 
                      as.character(gr.snv.enoc$ref), 
                      as.character(gr.snv.hgsc$ref))
  gr.snv.all$ALT <- c(as.character(unlist(gr.snv.pcawg$ALT)), 
                      as.character(gr.snv.ccoc$alt), 
                      as.character(gr.snv.enoc$alt), 
                      as.character(gr.snv.hgsc$alt))
  if(filterCoding){
    .cds <- cds
    if(filterSpliceSites){
      start(.cds) <- start(.cds)-5
      end(.cds) <- end(.cds)+5
    }
    gr.snv.all <- gr.snv.all[countOverlaps(query = gr.snv.all, subject = .cds) == 0]
  }
  
  seqlevels(gr.snv.all) <- paste0('chr',c(1:22,'X'))
  
  return(gr.snv.all)
}
