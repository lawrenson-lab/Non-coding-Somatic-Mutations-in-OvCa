.ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = 'grch37.ensembl.org')
.ensembl <- useDataset('hsapiens_gene_ensembl',mart=.ensembl)
.filters <- listFilters(.ensembl)
.attributes <- listAttributes(.ensembl)
bm <- getBM(attributes=c('ensembl_gene_id', 'description', 'gene_biotype', 'hgnc_symbol', 
                         'entrezgene_id', 'chromosome_name', 'start_position', 'end_position', 
                         'strand'),
            filters = 'chromosome_name',
            values = c(1:22,'X'),
            mart = .ensembl)

bm.genes <- GRanges(seqnames = paste0('chr',bm$chromosome_name), 
                    ranges = IRanges(start = bm$start_position, end = bm$end_position), 
                    strand = bm$strand, 
                    ensembl_gene_id = bm$ensembl_gene_id,
                    description = bm$description,
                    gene_biotype = bm$gene_biotype,
                    hgnc_symbol = bm$hgnc_symbol,
                    entrezgene = bm$entrezgene)

bm.tss <- GRanges(seqnames = paste0('chr',bm$chromosome_name), 
                  ranges = IRanges(start = bm$start_position, width = 1), 
                  strand = bm$strand, 
                  ensembl_gene_id = bm$ensembl_gene_id,
                  description = bm$description,
                  gene_biotype = bm$gene_biotype,
                  hgnc_symbol = bm$hgnc_symbol,
                  entrezgene = bm$entrezgene,
                  start_position = bm$start_position,
                  end_position = bm$end_position)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb, columns=c("gene_id"))
promoters <- promoters(txdb, upstream = 1000, downstream = 100, columns=c("gene_id"))
introns <- unlist(intronsByTranscript(txdb))
mcols(introns)$gene_id <- names(introns)
exons <- exons(txdb, columns = c('exon_id','gene_id'))
exons <- exons[which(unlist(lapply(X = mcols(exons)$gene_id, FUN = length))>0)]
cds <- cds(txdb)




getLocation <- function(gr){
  gr.tmp <- gr
  n <- length(gr)
  x <- c(reduce(genes, ignore.strand = T), reduce(promoters, ignore.strand = T))
  gr <- split(gr, countOverlaps(gr, x)==0)
  gr.intergenic <- gr$'TRUE'
  gr <- gr$'FALSE'
  gr <- split(gr, countOverlaps(gr, reduce(promoters, ignore.strand=T))==0)
  gr.promoter <- gr$'FALSE'
  gr <- gr$'TRUE'
  gr <- split(gr, countOverlaps(gr, reduce(introns, ignore.strand=T))==0)
  gr.intronic <- gr$'FALSE'
  gr.exonic <- gr$'TRUE'
  
  location <- rep(NA, length(gr.tmp))
  location[countOverlaps(gr.tmp, gr.intergenic)>0] <- 'intergenic'
  location[countOverlaps(gr.tmp, gr.promoter)>0] <- 'promoter'
  location[countOverlaps(gr.tmp, gr.intronic)>0] <- 'intronic'
  if(length(gr.exonic)>0)
    location[countOverlaps(gr.tmp, gr.exonic)>0] <- 'exonic'
  return(location)
}

getRawRNASeqFeatureCounts <- function(proteinCodingOnly = TRUE){
  x <- read.table('data/typecounts2.txt', header = T)
  tissue <- c(
    rep('CCOC',5),
    rep('EnOC',4),
    rep('HGSC',5),
    rep('MOC',5),
    rep('FTSEC',5)
  )
  patientId <- c(
    c('0241','0511','3172','3547','3588'),
    c('0291','0381','0630','0703'),
    c('0229','0270','0429','0550','0561'),
    c('0230','0380','0464','0652','3724'),
    c('1902','1907','1917','1939','1961')
  )
  colnames(x)[7:30] <- paste(tissue, patientId, sep = '_')
  rownames(x) <- x$Geneid
  
  if(proteinCodingOnly){
    x <- x[substr(x$Geneid, 0, 15) %in% bm$ensembl_gene_id[bm$gene_biotype=='protein_coding'],]
  }
  return(x)
}

read.dbaContrastRegions <- function(){
  x.ccoc <- read.table('Box/tissueChIP/results/dbaContrast_ccoc.tsv', header=T)
  gr.ccoc <- GenomicRanges::GRanges(seqnames = x.ccoc$seqnames, 
                                    ranges = IRanges::IRanges(start = x.ccoc$start, end = x.ccoc$end),
                                    strand = x.ccoc$strand, 
                                    Conc = x.ccoc$Conc, Conc_TRUE = x.ccoc$Conc_CCOC, Conc_FALSE = x.ccoc$Conc_NotCCOC, 
                                    Fold = x.ccoc$Fold, p.value = x.ccoc$p.value, FDR = x.ccoc$FDR)
  
  x.enoc <- read.table('Box/tissueChIP/results/dbaContrast_enoc.tsv', header=T)
  gr.enoc <- GenomicRanges::GRanges(seqnames = x.enoc$seqnames, 
                                    ranges = IRanges::IRanges(start = x.enoc$start, end = x.enoc$end),
                                    strand = x.enoc$strand, 
                                    Conc = x.enoc$Conc, Conc_TRUE = x.enoc$Conc_EnOC, Conc_FALSE = x.enoc$Conc_NotEnOC, 
                                    Fold = -x.enoc$Fold, p.value = x.enoc$p.value, FDR = x.enoc$FDR)
  
  x.hgsc <- read.table('Box/tissueChIP/results/dbaContrast_hgsc.tsv', header=T)
  gr.hgsc <- GenomicRanges::GRanges(seqnames = x.hgsc$seqnames, 
                                    ranges = IRanges::IRanges(start = x.hgsc$start, end = x.hgsc$end),
                                    strand = x.hgsc$strand, 
                                    Conc = x.hgsc$Conc, Conc_TRUE = x.hgsc$Conc_HGSOC, Conc_FALSE = x.hgsc$Conc_NotHGSC, 
                                    Fold = -x.hgsc$Fold, p.value = x.hgsc$p.value, FDR = x.hgsc$FDR)
  
  x.muoc <- read.table('Box/tissueChIP/results/dbaContrast_muoc.tsv', header=T)
  gr.muoc <- GenomicRanges::GRanges(seqnames = x.muoc$seqnames, 
                                    ranges = IRanges::IRanges(start = x.muoc$start, end = x.muoc$end),
                                    strand = x.muoc$strand, 
                                    Conc = x.muoc$Conc, Conc_TRUE = x.muoc$Conc_MOC, Conc_FALSE = x.muoc$Conc_NotMOC, 
                                    Fold = -x.muoc$Fold, p.value = x.muoc$p.value, FDR = x.muoc$FDR)
  
  gr2.ccoc <- gr.ccoc[countOverlaps(gr.ccoc, c(gr.enoc, gr.hgsc, gr.muoc))==0]
  gr2.enoc <- gr.enoc[countOverlaps(gr.enoc, c(gr.ccoc, gr.hgsc, gr.muoc))==0]
  gr2.hgsc <- gr.hgsc[countOverlaps(gr.hgsc, c(gr.ccoc, gr.enoc, gr.muoc))==0]
  gr2.muoc <- gr.muoc[countOverlaps(gr.muoc, c(gr.ccoc, gr.enoc, gr.hgsc))==0]
  
  mcols(gr2.ccoc)$histotype <- 'CCOC'
  mcols(gr2.enoc)$histotype <- 'EnOC'
  mcols(gr2.hgsc)$histotype <- 'HGSC'
  mcols(gr2.muoc)$histotype <- 'MuOC'
  
  gr <- c(gr2.ccoc, gr2.enoc, gr2.hgsc, gr2.muoc)
  return(gr)
}

read.histotypeSpecificRegions <- function(c1, c2){
  gr.ccoc <- import.bedGraph('results/specificToCCOC_all.bedGraph')
  gr.enoc <- import.bedGraph('results/specificToEnOC_all.bedGraph')
  gr.hgsc <- import.bedGraph('results/specificToHGSOC_all.bedGraph')
  gr.muoc <- import.bedGraph('results/specificToMOC_all.bedGraph')
  gr.ccoc <- reduce(gr.ccoc[mcols(gr.ccoc)$score >= c1])
  gr.enoc <- reduce(gr.enoc[mcols(gr.enoc)$score >= c1])
  gr.hgsc <- reduce(gr.hgsc[mcols(gr.hgsc)$score >= c1])
  gr.muoc <- reduce(gr.muoc[mcols(gr.muoc)$score >= c1])
  mcols(gr.ccoc)$histotype <- 'CCOC'
  mcols(gr.enoc)$histotype <- 'EnOC'
  mcols(gr.hgsc)$histotype <- 'HGSC'
  mcols(gr.muoc)$histotype <- 'MuOC'
  mcols(gr.ccoc)$foldPositive <- T
  mcols(gr.enoc)$foldPositive <- T
  mcols(gr.hgsc)$foldPositive <- T
  mcols(gr.muoc)$foldPositive <- T
  gr <- c(gr.ccoc, gr.enoc, gr.hgsc, gr.muoc)
  
  gr.ccoc <- import.bedGraph('results/specificallyMissingInCCOC_all.bedGraph')
  gr.enoc <- import.bedGraph('results/specificallyMissingInEnOC_all.bedGraph')
  gr.hgsc <- import.bedGraph('results/specificallyMissingInHGSOC_all.bedGraph')
  gr.muoc <- import.bedGraph('results/specificallyMissingInMOC_all.bedGraph')
  gr.ccoc <- reduce(gr.ccoc[mcols(gr.ccoc)$score >= c2])
  gr.enoc <- reduce(gr.enoc[mcols(gr.enoc)$score >= c2])
  gr.hgsc <- reduce(gr.hgsc[mcols(gr.hgsc)$score >= c2])
  gr.muoc <- reduce(gr.muoc[mcols(gr.muoc)$score >= c2])
  mcols(gr.ccoc)$histotype <- 'CCOC'
  mcols(gr.enoc)$histotype <- 'EnOC'
  mcols(gr.hgsc)$histotype <- 'HGSC'
  mcols(gr.muoc)$histotype <- 'MuOC'
  mcols(gr.ccoc)$foldPositive <- F
  mcols(gr.enoc)$foldPositive <- F
  mcols(gr.hgsc)$foldPositive <- F
  mcols(gr.muoc)$foldPositive <- F
  gr <- c(gr, gr.ccoc, gr.enoc, gr.hgsc, gr.muoc)
  return(gr)
}

getTADs <- function(extend = FALSE){
  hg19.len <- read.table(file = 'data/hg19.len')
  gr.tads <- import.bed(con = 'data/hesc1_TAD_domains_hg19.bed')
  gr.tads <- sort(gr.tads)
  if(extend){
    grl.tads <- split(x = gr.tads, f = seqnames(gr.tads))
    x <- foreach(i=1:length(grl.tads), .combine = c)%do%{
      .gr.tads <- grl.tads[[i]]
      y <- foreach(j=1:length(.gr.tads), .combine = c)%do%{
        z <- .gr.tads[j]
        if(j == 1){
          start(z) <- 1
        }else if(j == length(.gr.tads)){
          end(z) <- hg19.len$V2[as.character(hg19.len$V1) == as.character(seqnames(z))]
        }else{
          start(z) <- end(.gr.tads[j-1])
          end(z) <- start(.gr.tads[j+1])
        }
        return(z)
      }
      return(y)
    }
    gr.tads <- x
  }
  return(gr.tads)
}

getPromoters <- function(mart){
  bm.promoters <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'chromosome_name', 'transcript_start', 'strand'), mart = mart)
  bm.promoters$chromosome_name <- paste0('chr', bm.promoters$chromosome_name)
  colnames(bm.promoters)[3:4] <- c('seqnames', 'start')
  bm.promoters$end <- ifelse(test = bm.promoters$strand == 1, yes = bm.promoters$start + 100, no = bm.promoters$start + 1000)
  bm.promoters$start <- ifelse(test = bm.promoters$strand == 1, yes = bm.promoters$start - 1000, no = bm.promoters$start - 100)
  bm.promoters$strand <- ifelse(test = bm.promoters$strand == 1, yes = '+', no = '-')
  gr.promoters <- GRanges(bm.promoters)
  return(gr.promoters)
}

getOCMapping <- function(p.adjust.cutoff = 0.05, rPrime.cutoff = 0.5, mart){
  oc.seqnames <- paste0('chr',c(1:22,'X'))
  mapping <- read.table(file = 'results/mapEnahcner2Gene_genomeWide_TADs.txt', header = T, sep = '\t', quote = "")
  mapping <- mapping[mapping$p.adjust <= p.adjust.cutoff & mapping$r > 0 & mapping$rPrime > rPrime.cutoff,]
  seqnames <- unlist(strsplit(x = as.character(mapping$enhancer), split = ':'))
  seqnames <- seqnames[seq(from = 1, to = length(seqnames), by = 2)]
  range <- unlist(strsplit(x = as.character(mapping$enhancer), split = ':'))
  range <- range[seq(from = 2, to = length(range), by = 2)]
  range <- unlist(strsplit(x = range, split = '-'))
  start <- as.numeric(range[seq(from = 1, to = length(range), by = 2)])
  end <- as.numeric(range[seq(from = 2, to = length(range), by = 2)])
  mapping <- data.frame(seqnames, start, end, ensembl_gene_id = mapping$ensembl_gene_id)
  gr.mapping <- GRanges(mapping)
  gr.promoters <- getPromoters(mart = mart)
  mcols(gr.promoters) <- data.frame(ensembl_gene_id = gr.promoters$ensembl_gene_id)
  gr.mapping <- c(gr.mapping, gr.promoters)
  gr.mapping <- gr.mapping[seqnames(gr.mapping) %in% oc.seqnames]
  seqlevels(gr.mapping) <- oc.seqnames
  gr.mapping <- sort(gr.mapping)
  return(gr.mapping)
}
