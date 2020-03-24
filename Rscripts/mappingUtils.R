library(biomaRt)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(edgeR)
library(doMC)
library(foreach)


source('Rscripts/utils.R')
source('Rscripts/layer_IO.R')
source('Rscripts/layer_business.R')

mapping.init <- function(){
  load('RData/diffBind_session03_step02.RData')
  rm(list = 'mask')
  dba.DB <<- dba.peakset(DBA = dba, peaks = 1:20, minOverlap = 1, bRetrieve = T)
  x <- readRawCounts()
  x <- filterGenes(x)
  x <- filterSamples(x)
  normalizedCounts <<- getNormalizedCounts(x[,7:ncol(x)], x$Length, F)
  
  bm.genes <<- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'chromosome_name',
                     values = c(1:22,'X'),
                     mart = .ensembl)
  
  bm.tss <<- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'transcription_start_site'),
                   filters = 'chromosome_name',
                   values = c(1:22,'X'),
                   mart = .ensembl)
  
  gr.tads <<- getTADs(extend = TRUE)
}

mapEnhancerToGene <- function(gr.enhancer, ensembl_gene_id, type = 'within'){
  gr.enhancer <- subsetByOverlaps(dba.DB, gr.enhancer, type = type)
  row <- substr(x = rownames(normalizedCounts), start = 0, stop = 15) %in% ensembl_gene_id
  
  chipSeq <- mcols(gr.enhancer)
  chipSeq <- chipSeq[,-8]
  df <- data.frame(x = as.matrix(chipSeq)[1,], y = normalizedCounts[row,], 
                   histotype = c(rep('CCOC',5),rep('EnOC',4),rep('HGSOC',5),rep('MOC',5)))
  return(df)
}

mapEnhancerToGenes <- function(gr.enhancer, type = 'within', cor.method = 'spearman'){
  gr.enhancer <- subsetByOverlaps(query = dba.DB, subject = gr.enhancer, type = type)
  .gr.tads <- subsetByOverlaps(query = gr.tads, subject = gr.enhancer)
  ensembl_gene_ids <- unique(bm.tss$ensembl_gene_id[paste0('chr',bm.tss$chromosome_name) == as.character(seqnames(.gr.tads)) &
                                                      min(start(.gr.tads)) <= bm.tss$transcription_start_site &
                                                      max(end(.gr.tads)) >= bm.tss$transcription_start_site])
  ensembl_gene_ids <- ensembl_gene_ids[ensembl_gene_ids %in% substr(x = rownames(normalizedCounts), start = 0, stop = 15)]
  if(length(ensembl_gene_ids) > 0){
    c <- foreach(i=1:length(ensembl_gene_ids), .combine = rbind)%do%{
      ensembl_gene_id <- ensembl_gene_ids[i]
      d <- min(abs(c(start(gr.enhancer) - bm.tss$transcription_start_site[bm.tss$ensembl_gene_id == ensembl_gene_id],
                     end(gr.enhancer) - bm.tss$transcription_start_site[bm.tss$ensembl_gene_id == ensembl_gene_id])))
      df <- mapEnhancerToGene(gr.enhancer, ensembl_gene_id, type = type)
      cor <- cor.test(df$x, df$y, method = cor.method, alternative = 'greater')
      r <- cor$estimate
      p.value <- cor$p.value
      df <- data.frame(seqnames = seqnames(gr.enhancer), start = start(gr.enhancer), end = end(gr.enhancer), ensembl_gene_id, d, cor.method, r, p.value)
      return(df)
    }
    df <- as.data.frame(c)
  }else{
    df <- NA
  }
  return(df)
}

mapGeneToEnhancers <- function(ensembl_gene_id, cor.method = 'spearman'){
  bm.tss.index <- which(bm.tss$ensembl_gene_id == ensembl_gene_id)
  .gr.tads <- gr.tads[as.character(seqnames(gr.tads)) %in% paste0('chr', bm.tss$chromosome_name[bm.tss.index]) &
                        start(gr.tads) <= min(bm.tss$transcription_start_site[bm.tss.index]) &
                        end(gr.tads) >= max(bm.tss$transcription_start_site[bm.tss.index])]
  gr.enhancers <- subsetByOverlaps(query = dba.DB, subject = .gr.tads, type = 'within')
  if(length(gr.enhancers) > 0){
    c <- foreach(i=1:length(gr.enhancers), .combine = rbind)%do%{
      gr.enhancer <- gr.enhancers[i]
      d <- min(abs(c(start(gr.enhancer) - bm.tss$transcription_start_site[bm.tss.index],
                     end(gr.enhancer) - bm.tss$transcription_start_site[bm.tss.index])))
      df <- mapEnhancerToGene(gr.enhancer = gr.enhancer, ensembl_gene_id = ensembl_gene_id, type = 'equal')
      cor <- cor.test(df$x, df$y, method = cor.method, alternative = 'greater')
      r <- cor$estimate
      p.value <- cor$p.value
      df <- data.frame(seqnames = seqnames(gr.enhancer), start = start(gr.enhancer), end = end(gr.enhancer), ensembl_gene_id, d, cor.method, r, p.value)
      return(df)
    }
    df <- as.data.frame(c)
  }else{
    df <- NA
  }
  return(df)
}

mapEnhancersToGenes <- function(.gr.tads, cor.method = 'spearman'){
  gr.enhancers <- subsetByOverlaps(query = dba.DB, subject = .gr.tads, type = 'within')
  ensembl_gene_ids <- unique(bm.tss$ensembl_gene_id[paste0('chr',bm.tss$chromosome_name) == as.character(seqnames(.gr.tads)) &
                                                      min(start(.gr.tads)) <= bm.tss$transcription_start_site &
                                                      max(end(.gr.tads)) >= bm.tss$transcription_start_site])
  ensembl_gene_ids <- ensembl_gene_ids[ensembl_gene_ids %in% substr(x = rownames(normalizedCounts), start = 0, stop = 15)]
  if(length(gr.enhancers)>0 & length(ensembl_gene_ids)>0){
    df <- foreach(i=1:length(gr.enhancers), .combine = rbind)%do%{
      gr.enhancer <- gr.enhancers[i]
      df <- foreach(ensembl_gene_id=ensembl_gene_ids, .combine = rbind)%do%{
        d <- min(abs(c(start(gr.enhancer) - bm.tss$transcription_start_site[bm.tss$ensembl_gene_id == ensembl_gene_id],
                       end(gr.enhancer) - bm.tss$transcription_start_site[bm.tss$ensembl_gene_id == ensembl_gene_id])))
        df <- mapEnhancerToGene(gr.enhancer = gr.enhancer, ensembl_gene_id = ensembl_gene_id, type = 'equal')
        cor <- cor.test(df$x, df$y, method = cor.method, alternative = 'greater')
        r <- cor$estimate
        p.value <- cor$p.value
        df <- data.frame(seqnames = seqnames(gr.enhancer), start = start(gr.enhancer), end = end(gr.enhancer), ensembl_gene_id, d, cor.method, r, p.value)
        return(df)
      }
      return(df)
    }
    return(df)
  }
  return(NA)
}




.plotEnhancerToGene <- function(data, enhancer.name, gene.name){
  cor.pearson <- cor.test(x = data$x, y = data$y, method = 'pearson', alternative = 'greater')
  cor.kendall <- cor.test(x = data$x, y = data$y, method = 'kendall', alternative = 'greater')
  cor.spearman <- cor.test(x = data$x, y = data$y, method = 'spearman', alternative = 'greater')
  p <- ggplot(data = data, mapping = aes(x, y, col=histotype)) + 
    geom_point() + 
    ggtitle(label = hgnc_symbol, subtitle = paste0("Pearson's correlation, r = ",round(cor.pearson$estimate,2),
                                                   ", p = ",round(cor.pearson$p.value,4),
                                                   "\nSpearman's correlation, r = ",round(cor.spearman$estimate,2),
                                                   ", p = ",round(cor.spearman$p.value,4),
                                                   "\nKendall's correlation, r = ",round(cor.kendall$estimate,2),
                                                   ", p = ",round(cor.kendall$p.value,4))) + 
    scale_color_brewer(palette = 'Set1') + 
    xlab(paste0(enhancer.name,'\nNormalized H3K27ac ChIP-Seq signal')) + 
    ylab(paste0(gene.name,' expression\nRNA-Seq (CPM)'))
  cor.pearson <- cor.test(x = log2(data$x), y = log2(data$y), method = 'pearson', alternative = 'greater')
  cor.kendall <- cor.test(x = log2(data$x), y = log2(data$y), method = 'kendall', alternative = 'greater')
  cor.spearman <- cor.test(x = log2(data$x), y = log2(data$y), method = 'spearman', alternative = 'greater')
  q <- ggplot(data = data, mapping = aes(x, y, col=histotype)) + 
    geom_point() + 
    ggtitle(label = hgnc_symbol, subtitle = paste0("Pearson's correlation, r = ",round(cor.pearson$estimate,2),
                                                   ", p = ",round(cor.pearson$p.value,4),
                                                   "\nSpearman's correlation, r = ",round(cor.spearman$estimate,2),
                                                   ", p = ",round(cor.spearman$p.value,4),
                                                   "\nKendall's correlation, r = ",round(cor.kendall$estimate,2),
                                                   ", p = ",round(cor.kendall$p.value,4))) + 
    scale_color_brewer(palette = 'Set1') + 
    xlab(paste0(enhancer.name,'\nNormalized H3K27ac ChIP-Seq signal')) + 
    ylab(paste0(gene.name,' expression\nRNA-Seq (CPM)')) + 
    scale_x_continuous(trans = 'log2') + 
    scale_y_continuous(trans = 'log2')
  grid.arrange(p, q, nrow = 1)
}

.plotGeneToEnhancers <- function(ensembl_gene_id){
  bm.tss.index <- which(bm.tss$ensembl_gene_id == ensembl_gene_id)
  .gr.tads <- gr.tads[as.character(seqnames(gr.tads)) %in% paste0('chr', bm.tss$chromosome_name[bm.tss.index]) &
                        start(gr.tads) <= min(bm.tss$transcription_start_site[bm.tss.index]) &
                        end(gr.tads) >= max(bm.tss$transcription_start_site[bm.tss.index])]
  gr.enhancers <- subsetByOverlaps(query = dba.DB, subject = .gr.tads, type = 'within')
  for(i in 1:length(gr.enhancers)){
    gr.enhancer <- gr.enhancers[i]
    enhancer.name <- paste0(seqnames(gr.enhancer),':',start(gr.enhancer),'-',end(gr.enhancer))
    data <- mapEnhancerToGene(gr.enhancer = gr.enhancer, ensembl_gene_id = ensembl_gene_id)
    gene.name <- hgnc_symbol
    enhancer.name <- paste0(seqnames(gr.enhancer),':',start(gr.enhancer),'-',end(gr.enhancer))
    print(.plotEnhancerToGene(data = data, enhancer.name = enhancer.name, gene.name = gene.name))
  }
}


mapTADbyGene <- function(hgnc_symbol, cor.method){
  ensembl_gene_id <- bm.genes$ensembl_gene_id[bm.genes$hgnc_symbol == hgnc_symbol]
  bm.tss.index <- which(bm.tss$ensembl_gene_id == ensembl_gene_id)
  .gr.tads <- gr.tads[as.character(seqnames(gr.tads)) %in% paste0('chr', bm.tss$chromosome_name[bm.tss.index]) &
                        start(gr.tads) <= min(bm.tss$transcription_start_site[bm.tss.index]) &
                        end(gr.tads) >= max(bm.tss$transcription_start_site[bm.tss.index])]
  ensembl_gene_ids <- unique(bm.tss$ensembl_gene_id[paste0('chr',bm.tss$chromosome_name) == as.character(seqnames(.gr.tads)) &
                                                      min(start(.gr.tads)) <= bm.tss$transcription_start_site &
                                                      max(end(.gr.tads)) >= bm.tss$transcription_start_site])
  ensembl_gene_ids <- ensembl_gene_ids[ensembl_gene_ids %in% substr(x = rownames(normalizedCounts), start = 0, stop = 15)]
  for(ensembl_gene_id in ensembl_gene_ids){
    print(ensembl_gene_id)
    df.enhancerGeneMap <- mapGeneToEnhancers(ensembl_gene_id = ensembl_gene_id, cor.method = cor.method)
    write.table(x = df.enhancerGeneMap, file = paste0('results/enhnacerGeneMappingDB.txt'), append = T, quote = F, sep = '\t', row.names = F, col.names = F)
    is.first <- F
  }
}
