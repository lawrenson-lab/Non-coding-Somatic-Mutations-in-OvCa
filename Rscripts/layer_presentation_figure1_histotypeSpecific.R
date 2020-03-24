library(ggplot2)
library(scales)
library(foreach)
library(doMC)

getIndex <- function(gr){
  histotypes <- unique(gr$histotype)
  x <- foreach(j=1:length(histotypes), .combine = rbind)%dopar%{
    .gr <- gr[gr$histotype == histotypes[j]]
    gr.pos <- .gr[.gr$Fold > 0]
    gr.neg <- .gr[.gr$Fold < 0]
    gr.pos <- gr.pos[order(-mcols(gr.pos)[,2])]
    gr.neg <- gr.neg[order(-mcols(gr.neg)[,3])]
    .gr <- c(gr.pos, gr.neg)
    index <- findOverlaps(query = .gr, subject = dba.DB, type = 'equal', select = 'first')
    group <- c(rep('Positive', length(gr.pos)), rep('Negative', length(gr.neg)))
    histotype <- histotypes[j]
    index <- data.frame(index, group, histotype)
    return(index)
  }
  return(x)
}

plotLocation <- function(gr, legend.position = 'none'){
  gr$location <- getLocation(gr)
  df <- as.data.frame(gr)
  df$n <- 1
  df <- ddply(df, .(group, location), summarize, n = sum(n))
  print(df)
  df$location <- factor(df$location)
  p <- ggplot(data = df, mapping = aes(x = group, y = n, col=location, fill=location))
  p + geom_bar(stat = 'identity', position='dodge') + coord_flip() +
    scale_y_continuous(labels = comma, trans = 'reverse') + expand_limits(x = 0) + 
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          legend.position = legend.position) 
}