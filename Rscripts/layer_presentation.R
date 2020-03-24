w <- 89/2.54/10

plotExponentialDist <- function(ex, p.value = 0.05, main = '', xlab = ''){
  fit <- fitdistr(x = ex, densfun = 'exponential')
  hist(ex, freq = FALSE, breaks = 100, xlim = c(0, quantile(ex, 0.99)), main = main, xlab = xlab)
  curve(dexp(x, rate = fit$estimate), col = "red", add = TRUE, lwd = 2)
  abline(v = qexp(1-p.value, rate = fit$estimate), lty = 2, col = 'red')
}

plotExponentialDist2 <- function(ex, p.value = 0.05, main = ''){
  fit <- fitdistr(x = ex, densfun = 'exponential')
  hist(ex, freq = FALSE, main = main)
  curve(dexp(x, rate = fit$estimate), col = "red", add = TRUE, lwd = 2)
  abline(v = qexp(1-p.value, rate = fit$estimate), lty = 2, col = 'red')
}

plotNormalDist <- function(ex, p.value = 0.05, main = ''){
  fit <- fitdistr(x = ex, densfun = 'normal')
  hist(ex, freq = FALSE, breaks = 100, xlim = c(quantile(ex, 0), quantile(ex, 0.99)), main = main)
  curve(dnorm(x, mean = fit$estimate[1], sd = fit$estimate[2]), col = "red", add = TRUE, lwd = 2)
  abline(v = qnorm(1-p.value, mean = fit$estimate[1], sd = fit$estimate[2]), lty = 2, col = 'red')
}

plotLogisticDist <- function(ex, p.value = 0.05, main = ''){
  fit <- fitdistr(x = ex, densfun = 'logistic')
  hist(ex, freq = FALSE, breaks = 200, xlim = c(0, quantile(ex, 0.99)), main = main)
  curve(dlogis(x, location = fit$estimate[1], scale = fit$estimate[2]), col = "red", add = TRUE, lwd = 2)
  abline(v = qlogis(1-p.value, location = fit$estimate[1], scale = fit$estimate[2]), lty = 2, col = 'red')
}

plotPoissonDist <- function(ex, p.value = 0.05, main = ''){
  fit <- fitdistr(x = ex, densfun = 'poisson')
  hist(ex, freq = FALSE, breaks = 200, xlim = c(0, quantile(ex, 0.99)), main = main)
  points(dpois(1:4000, lambda = fit$estimate), col = "red", lwd = 2)
  abline(v = qpois(1-p.value, lambda = fit$estimate), lty = 2, col = 'red')
}

plotCauchyDist <- function(ex, p.value = 0.05, main = ''){
  fit <- fitdistr(x = ex, densfun = 'cauchy')
  hist(ex, freq = FALSE, breaks = 200, xlim = c(0, quantile(ex, 0.99)), main = main)
  curve(dcauchy(x, location = fit$estimate[1], scale = fit$estimate[2]), col = "red", lwd = 2, add = T)
  abline(v = qcauchy(1-p.value, location = fit$estimate[1], scale = fit$estimate[2]), lty = 2, col = 'red')
}

plotDips <- function(.gr.narrowPeaks, gr.bigWig, order = 250, plot = F, main = ''){
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
  if(plot){
    df <- as.data.frame(z)
    colnames(df) <- c('x','y')
    df$y2 <- y2
    p <- ggplot(data = df, mapping = aes(x = x, y = y))
    print(p + xlab(seqnames(.gr.narrowPeaks)) + ylab('Fold change') + ggtitle(label = main, subtitle = .name) +
            geom_line(col = col[1]) +
            geom_line(mapping = aes(x = x, y = y2), col = col[2], size = 1) +
            geom_vline(xintercept = z[iLocalMin,1], col = col[3], size = 1, linetype = 'dotted') +
            geom_vline(xintercept = z[iLocalMax,1], col = col[4], size = 1, linetype = 'dotted'))
  }

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

plotManhattan <- function(gr.H3K27ac, pvalues, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08)){
  par(las = 2)
  seqnames <- paste0('chr',c(1:22,'X'))
  seqlevels(gr.H3K27ac) <- seqnames
  df <- data.frame(SNP = paste(seqnames(gr.H3K27ac),start(gr.H3K27ac),sep = '_'), 
                   CHR = as.numeric(seqnames(gr.H3K27ac)), BP = start(gr.H3K27ac), 
                   P = pvalues)
  manhattan(x = df, col = c('#DE47AB', '#72BE97', '#F7F797', '#7C749B', '#E85726',
                            '#B395F8', '#DC8747', '#96D53D', '#DC85EE', '#7D32B3',
                            '#88DB68', '#78AAF1', '#D9C6CA', '#336C80', '#F7CA44',
                            '#32C7C7', '#D4C5F2', '#995493', '#F88B78', '#475ECC',
                            '#E0BD8C', '#9E2800', '#F2BBD2'), chrlabs = seqnames, 
            suggestiveline = suggestiveline, genomewideline = genomewideline)
  abline(h = .x1+.x2, col='black', lty = 2, lwd = 2)
}

plot.chipSeqSignal <- function(gr, score = DBA_SCORE_TMM_MINUS_FULL_CPM){
  if(is.na(dba.DB)){
    load('RData/diffBind_session03_step02.RData')
    dba <- dba.count(DBA = dba, peaks = NULL, score = score)
    dba.DB <<- dba.peakset(DBA = dba, peaks = 1:20, minOverlap = 1, bRetrieve = T)
  }
  gr <- subsetByOverlaps(query = dba.DB, subject = gr, type = 'within')
  regionName <- paste0(as.character(seqnames(gr)), ":", start(gr), '-', end(gr))
  chipSeqSignal <- t(as.matrix(mcols(gr)))
  histotype <- substr(x = rownames(chipSeqSignal), start = 1, stop = nchar(rownames(chipSeqSignal))-2)
  df <- data.frame(histotype, chipSeqSignal)
  p <- ggplot(data = df, mapping = aes(x = histotype, y = chipSeqSignal, col = histotype))
  p + geom_boxplot(show.legend = F, outlier.colour = NA) + scale_color_brewer(palette = 'Set1') + ylab(paste0(regionName,'\nH3K27ac ChIP-Seq signal (CPM)')) + theme(axis.title.x = element_blank()) + 
    expand_limits(x = 0, y = 0) + geom_jitter(width = 0.1, show.legend = F)
}

plot.rnaSeq <- function(hgnc_symbol){
  if(is.na(ensembl_hgnc_mapping)){
    .ensembl.hg38 <- useMart("ensembl")
    .ensembl.hg38 <- useDataset('hsapiens_gene_ensembl',mart=.ensembl.hg38)
    ensembl_hgnc_mapping <<- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                                   filters = 'chromosome_name',
                                   values = c(1:22,'X'),
                                   mart = .ensembl.hg38)
  }
  if(is.na(normalizedCounts)){
    x <- readRawCounts()
    x <- filterGenes(x)
    x <- filterSamples(x)
    normalizedCounts <<- getNormalizedCounts(x[,7:ncol(x)], x$Length, F)
  }
  ylab <- paste(hgnc_symbol,' expression\nRNA-Seq (CPM)')
  row <- startsWith(rownames(normalizedCounts),ensembl_hgnc_mapping$ensembl_gene_id[ensembl_hgnc_mapping$hgnc_symbol == hgnc_symbol])
  histotype <- substr(x = colnames(normalizedCounts), start = 1, stop = nchar(colnames(normalizedCounts))-2)
  data <- data.frame(rnaSeq = normalizedCounts[row,], histotype)
  p <- ggplot(data = data, mapping = aes(col = histotype, y = rnaSeq, x = histotype))
  p + geom_boxplot(show.legend = F, outlier.colour = NA) + scale_color_brewer(palette = 'Set1') + ylab(ylab) + expand_limits(x = 0, y = 0) + geom_jitter(width = 0.1, show.legend = F) + 
    theme(axis.title.x = element_blank())
}

plot.chipSeq_rnaSeq_scatterplot <- function(gr, score = DBA_SCORE_TMM_MINUS_FULL_CPM, hgnc_symbol){
  if(is.na(dba.DB)){
    load('RData/diffBind_session03_step02.RData')
    dba <- dba.count(DBA = dba, peaks = NULL, score = score)
    dba.DB <<- dba.peakset(DBA = dba, peaks = 1:20, minOverlap = 1, bRetrieve = T)
  }
  if(is.na(ensembl_hgnc_mapping)){
    .ensembl.hg38 <- useMart("ensembl")
    .ensembl.hg38 <- useDataset('hsapiens_gene_ensembl',mart=.ensembl.hg38)
    ensembl_hgnc_mapping <<- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                                   filters = 'chromosome_name',
                                   values = c(1:22,'X'),
                                   mart = .ensembl.hg38)
  }
  if(is.na(normalizedCounts)){
    x <- readRawCounts()
    x <- filterGenes(x)
    x <- filterSamples(x)
    normalizedCounts <<- getNormalizedCounts(x[,7:ncol(x)], x$Length, F)
  }
  
  gr <- subsetByOverlaps(query = dba.DB, subject = gr, type = 'within')
  regionName <- paste0(as.character(seqnames(gr)), ":", start(gr), '-', end(gr))
  chipSeqSignal <- t(as.matrix(mcols(gr)))[-8]
  ylab <- paste(hgnc_symbol,' expression\nRNA-Seq (CPM)')
  row <- startsWith(rownames(normalizedCounts),ensembl_hgnc_mapping$ensembl_gene_id[ensembl_hgnc_mapping$hgnc_symbol == hgnc_symbol])
  histotype <- substr(x = colnames(normalizedCounts), start = 1, stop = nchar(colnames(normalizedCounts))-2)
  data <- data.frame(rnaSeq = normalizedCounts[row,], histotype, chipSeqSignal)
  p <- ggplot(data = data, mapping = aes(col = histotype, y = rnaSeq, x = chipSeqSignal))
  p + geom_point(show.legend = F) + scale_color_brewer(palette = 'Set1') + ylab(ylab) + expand_limits(x = 0, y = 0) + 
    xlab(paste0(regionName,'\nH3K27ac ChIP-Seq signal (CPM)')) + scale_x_log10(labels = comma) + scale_y_log10(labels = comma)
}
