readRawCounts <- function(){
  x <- read.table(file = 'data/typecounts2.txt', header = T)
  tissue <- c(
    rep('CCOC',5),
    rep('EnOC',4),
    rep('HGSOC',5),
    rep('MOC',5),
    rep('FTSEC',5)
  )
  patientId <- c(1:5,1:2,4:5,1:5,1:5,1:5)
  colnames(x)[7:30] <- paste(tissue, patientId, sep = '-')
  rownames(x) <- x$Geneid
  return(x)
}

readChIPSeqQC <- function(){
  x <- read.table(file = 'results/qc.txt', header = T, sep = '\t')
  x <- x[c(1,3:6,9,12,15:17,19:21,18,22:27),]
  x$Histotype <- c(rep('CCOC',5),
                   rep('EnOC',5),
                   rep('HGSOC',5),
                   rep('MOC',5))
  return(x)
}

pcawg.getAliquotIdsByTumorSubtype <- function(tumorSubtype){
  x <- read.table(file = 'data/PCAWG/whitelisted_tumour_aliqouts.20160831.txt', header = T, sep = '\t')
  if(tumorSubtype == 'Ovary-AdenoCA'){
    return(x$aliquot_id[x$dcc_project_code %in% c('OV-AU', 'OV-US')])
  }else if(tumorSubtype == 'Prost-AdenoCA'){
    return(x$aliquot_id[x$dcc_project_code %in% c('PRAD-CA', 'EOPC-DE', 'PRAD-UK', 'PRAD-US')])
  }
  return(NA)
}

pcawg.getAliquotsByTumorSubtype <- function(tumorSubtype){
  x <- read.table(file = 'data/PCAWG/whitelisted_tumour_aliqouts.20160831.txt', header = T, sep = '\t')
  if(tumorSubtype == 'Ovary-AdenoCA'){
    return(x[x$dcc_project_code %in% c('OV-AU', 'OV-US'),])
  }else if(tumorSubtype == 'Prost-AdenoCA'){
    return(x[x$dcc_project_code %in% c('PRAD-CA', 'EOPC-DE', 'PRAD-UK', 'PRAD-US'),])
  }
  return(NA)
}

pcawg.getAliquotsByHistologyAbbreviation <- function(histology_abbreviation){
  x <- read.table(file = 'data/PCAWG/whitelisted_tumour_aliqouts.20160831.txt', header = T, sep = '\t')
  y <- read.table(file = 'data/PCAWG/PCAWG-13_10 Molecular subtypes and Clinical correlates/pcawg_specimen_histology_August2016_v9.tsv', header = T, sep = '\t', quote = "", comment.char = "")
  aliqouts <- x[x$icgc_specimen_id %in% y$X..icgc_specimen_id[y$histology_abbreviation == histology_abbreviation],]
  return(aliqouts)
}

pcawg.getAliquots <- function(){
  x <- read.table(file = 'data/PCAWG/whitelisted_tumour_aliqouts.20160831.txt', header = T, sep = '\t')
  y <- read.table(file = 'data/PCAWG/PCAWG-13_10 Molecular subtypes and Clinical correlates/pcawg_specimen_histology_August2016_v9.tsv', header = T, sep = '\t', quote = "", comment.char = "")
  aliqouts <- merge(x = x, y = y, by.x = 'icgc_specimen_id', by.y = 'X..icgc_specimen_id')
  return(aliqouts)
}

pcawg.getHistologyAbbreviations <- function(){
  x <- read.table(file = 'data/PCAWG/PCAWG-13_10 Molecular subtypes and Clinical correlates/pcawg_donor_subtype_cohort_list.tsv', header = T, sep = '\t', nrows = 34, skip = 3)
  return(x$histology_abbreviation)
}

pcawg.getSNVByAliquotId <- function(aliquot_id){
  path <- 'data/PCAWG/final_consensus_12oct_passonly/snv_mnv/'
  file <- paste0(path, aliquot_id,'.consensus.20160830.somatic.snv_mnv.vcf.gz')
  vcf <- readVcf(file = file, genome = 'hg19')
  gr.tmp <- rowRanges(vcf)
  gr.tmp <- GRanges(seqnames = paste0('chr', seqnames(gr.tmp)), 
                    ranges = ranges(gr.tmp), 
                    strand = strand(gr.tmp),
                    paramRangeID = mcols(gr.tmp)$paramRangeID,
                    REF = mcols(gr.tmp)$REF,
                    ALT = mcols(gr.tmp)$ALT,
                    QUAL = mcols(gr.tmp)$QUAL,
                    FILTER = mcols(gr.tmp)$FILTER,
                    '1000genomes_AF' = info(vcf)[,which(names(info(vcf))=='1000genomes_AF')],
                    '1000genomes_ID' = info(vcf)[,which(names(info(vcf))=='1000genomes_ID')],
                    Callers = info(vcf)[,which(names(info(vcf))=='Callers')],
                    NumCallers = info(vcf)[,which(names(info(vcf))=='NumCallers')],
                    VAF = info(vcf)[,which(names(info(vcf))=='VAF')],
                    cosmic = info(vcf)[,which(names(info(vcf))=='cosmic')],
                    dbsnp = info(vcf)[,which(names(info(vcf))=='dbsnp')],
                    repeat_masker = info(vcf)[,which(names(info(vcf))=='repeat_masker')],
                    t_alt_count = info(vcf)[,which(names(info(vcf))=='t_alt_count')],
                    t_ref_count = info(vcf)[,which(names(info(vcf))=='t_ref_count')],
                    dbsnp_somatic = info(vcf)[,which(names(info(vcf))=='dbsnp_somatic')],
                    Variant_Classification = info(vcf)[,which(names(info(vcf))=='Variant_Classification')])
  return(gr.tmp)
}

pcawg.getInDelsByAliquotId <- function(aliquot_id){
  path <- 'data/PCAWG/final_consensus_12oct_passonly/indel/'
  file <- paste0(path, aliquot_id,'.consensus.20161006.somatic.indel.vcf.gz')
  vcf <- readVcf(file)
  gr.tmp <- rowRanges(vcf)
  gr.tmp <- GRanges(seqnames = paste0('chr', seqnames(gr.tmp)), 
                    ranges = ranges(gr.tmp), 
                    strand = strand(gr.tmp),
                    paramRangeID = mcols(gr.tmp)$paramRangeID,
                    REF = mcols(gr.tmp)$REF,
                    ALT = mcols(gr.tmp)$ALT,
                    QUAL = mcols(gr.tmp)$QUAL,
                    FILTER = mcols(gr.tmp)$FILTER,
                    '1000genomes_AF' = info(vcf)[,which(names(info(vcf))=='1000genomes_AF')],
                    '1000genomes_ID' = info(vcf)[,which(names(info(vcf))=='1000genomes_ID')],
                    Callers = info(vcf)[,which(names(info(vcf))=='Callers')],
                    NumCallers = info(vcf)[,which(names(info(vcf))=='NumCallers')],
                    VAF = info(vcf)[,which(names(info(vcf))=='VAF')],
                    cosmic = info(vcf)[,which(names(info(vcf))=='cosmic')],
                    dbsnp = info(vcf)[,which(names(info(vcf))=='dbsnp')],
                    repeat_masker = info(vcf)[,which(names(info(vcf))=='repeat_masker')],
                    t_alt_count = info(vcf)[,which(names(info(vcf))=='t_alt_count')],
                    t_ref_count = info(vcf)[,which(names(info(vcf))=='t_ref_count')],
                    dbsnp_somatic = info(vcf)[,which(names(info(vcf))=='dbsnp_somatic')],
                    Variant_Classification = info(vcf)[,which(names(info(vcf))=='Variant_Classification')])
  return(gr.tmp)
}


pcawg.getSNVsByHistologyAbbreviation <- function(histology_abbreviation = 'Ovary-AdenoCA'){
  aliquots <- pcawg.getAliquotsByHistologyAbbreviation(histology_abbreviation)
  gr <- foreach(aliquot_id=aliquots$aliquot_id, .combine = c) %do% {
    .gr1 <- pcawg.getSNVByAliquotId(aliquot_id = aliquot_id)
    .gr1$type <- 'SNV'
    .gr <- .gr1
    .gr$aliquot_id <- aliquot_id
    return(.gr)
  }
  return(gr)
}

pcawg.readMapEnhancerGene <- function(){
  x <- read.table(file = 'data/PCAWG/NOT_SET', header = F, sep = '\t', quote = "", as.is = T)
  target.genes <- x$V2
  x <- unlist(strsplit(x = x$V1, split = ':'))
  seqnames <- x[seq(1,length(x)-1, by = 2)]
  x <- x[seq(2,length(x), by = 2)]
  x <- unlist(strsplit(x = as.character(x), split = '-'))
  start <- as.numeric(x[seq(1,length(x)-1, by = 2)])
  end <- as.numeric(x[seq(2,length(x), by = 2)])
  gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end-1), target.genes = target.genes)
  return(gr)
}

pcawg.getSNVs <- function(histology_abbreviation, filterCoding = TRUE, filterSpliceSites = F){
  gr.snv.pcawg <- pcawg.getSNVsByHistologyAbbreviation(histology_abbreviation = histology_abbreviation)
  gr.snv.all <- granges(gr.snv.pcawg)
  gr.snv.all$id <- as.character(gr.snv.pcawg$aliquot_id)
  gr.snv.all$dataset <- rep('PCAWG',length(gr.snv.pcawg))
  gr.snv.all$REF <- as.character(gr.snv.pcawg$REF)
  gr.snv.all$ALT <- as.character(unlist(gr.snv.pcawg$ALT))
  if(filterCoding){
    .cds <- cds
    if(filterSpliceSites){
      start(.cds) <- start(.cds)-5
      end(.cds) <- end(.cds)+5
    }
    gr.snv.all <- gr.snv.all[countOverlaps(query = gr.snv.all, subject = .cds) == 0]
  }
  gr.snv.all <- gr.snv.all[seqlevels(gr.snv.all) %in% paste0('chr',c(1:22,'X'))]
  return(gr.snv.all)
}

shah.getCaseIdsByProject <- function(project = 'HGSC'){
  #reads Shah's data and returns the case ids for a specific project
  #available projects are: CCOC, ENOC, GCT and HGSC
  file <- paste0('data/OvCa.ShahLab.WGS/MUT_', project, '_.txt')
  x <- read.table(file = file, header = T, sep = '\t')
  return(unique(x$case_id))
}

shah.getSNVsByProject <- function(project = 'HGSC', remove.EnOC.outlier = F){
  #reads Shah's data and returns the mutations for a specific project
  #available projects are: CCOC, ENOC, GCT and HGSC
  file <- paste0('data/OvCa.ShahLab.WGS/MUT_', project, '_.txt')
  x <- read.table(file = file, header = T, sep = '\t')
  x <- x[nchar(as.character(x$ref)) == nchar(as.character(x$alt)),]
  gr <- GRanges(seqnames = paste0('chr', x$chromosome), ranges = IRanges(start = x$start, end = x$stop))
  mcols(gr) <- x[,!colnames(x) %in% c('chromosome','start','stop')]
  if(remove.EnOC.outlier){
    gr <- gr[!(gr$case_id == 'DG1285')]
  }
  return(gr)
}

writeMAFLITE <- function(gr, fileout, maxLines = NA){
  df <- as.data.frame(gr)
  df$ALT <- as.character(unlist(mcols(gr)$ALT))
  df <- data.frame(df$seqnames, df$start, df$end, df$REF, df$ALT)
  colnames(df) <- c('chr', 'start', 'end', 'ref_allele', 'alt_allele')
  df$ref_allele <- as.character(df$ref_allele)
  df$alt_allele <- as.character(df$alt_allele)
  if(is.na(maxLines) | nrow(df) <= maxLines){
    write.table(df, fileout, append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  }else{
    maxLines <- maxLines-1
    for(i in 0:floor((nrow(df)-1)/maxLines)){
      .df <- df[floor((1:nrow(df)-1)/maxLines) == i,]
      write.table(.df, paste0(substr(fileout, 1, nchar(fileout)-8),'_',i+1,'.maflite'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)
    }
  }
}

getBlacklist <- function(){
  blacklist <- import.bed(con = 'data/UCSCGenomeBrowser/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
  blacklist <- blacklist[seqnames(blacklist) %in% paste0('chr',c(1:22,'X'))]
  return(blacklist)
}

read.narrowPeak <- function(filename){
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
  gr <- import(con = filename, format = "BED", extraCols = extraCols_narrowPeak)
  return(gr)
}

read.IDROutput <- function(filename, numberOfRep = 2){
  gr <- NA
  if(numberOfRep == 2){
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", 
                              qValue = "numeric", summit = "integer", 
                              localIDR = "numeric", globalIDR = "numeric", 
                              rep1_chromStart = "integer", rep1_chromEnd = "integer", 
                              rep1_signalValue = "numeric", rep1_summit = "integer", 
                              rep2_chromStart = "integer", rep2_chromEnd = "integer", 
                              rep2_signalValue = "numeric", rep2_summit = "integer")
    gr <- import.bed(con = filename, extraCols = extraCols_narrowPeak)
  }
  return(gr)
}

import.bed.withCommentedLines <- function(file){
  .gr <- read.table(file = file, header = F, sep = '\t', quote = "", as.is = T, comment.char = '#', blank.lines.skip = T)
  .gr <- .gr[startsWith(x = .gr$V1, prefix = 'chr') | .gr$V1 %in% c(1:22,'X','Y'),]
  gr <- GRanges(seqnames = .gr$V1, ranges = IRanges(start = .gr$V2, end = .gr$V3))
  if(ncol(.gr)>3){mcols(gr)$name <- .gr$V4}
  if(ncol(.gr)>4){mcols(gr)$score <- .gr$V5}
  if(ncol(.gr)>5){mcols(gr)$strand <- .gr$V6}
  return(gr)
}

homer.readSuperEnhancers <- function(filename){
  df <- read.table(filename)
  colnames(df) <- c('name', 'seqnames', 'start', 'end', 'strand', 'normalizedTagCount', 'superEnahncerSlope', 
                    'peakScore', 'totalTags')#, 'controlTags', 'foldChange', 'pValue', 'clonalFoldChange')
  gr <- GRanges(seqnames = df$seqnames,
                ranges = IRanges(start = df$start, end = df$end),
                strand = '*', name = df$name,
                normalizedTagCount = df$normalizedTagCount, superEnahncerSlope = df$superEnahncerSlope,
                peakScore = df$peakScore, totalTags = df$totalTags)
  return(gr)
}

homer.readHistoneNFR <- function(filename){
  df <- read.table(filename)
  colnames(df) <- c('name', 'seqnames', 'start', 'end', 'strand', 'normalizedTagCount', 'regionSize', 'peakScore', 
                    'totalTags', 'controlTags', 'foldChange', 'pValue', 'clonalFoldChange')
  gr <- GRanges(seqnames = df$seqnames,
                ranges = IRanges(start = df$start, end = df$end),
                strand = '*', name = df$name,
                normalizedTagCount = df$normalizedTagCount, regionSize = df$regionSize, peakScore = df$peakScore, 
                totalTags = df$totalTags, controlTags = df$controlTags, foldChange = df$foldChange,
                pValue = df$pValue, clonalFoldChange = df$clonalFoldChange)
  return(gr)
}

read.somMutations <- function(filename){
  extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer", 
                 nSNV = "integer", mutRate = "numeric", localMutRate = "numeric", adjMutRate = "numeric", 
                 nSamples = "integer", dSamples = "numeric")
  gr <- import(con = filename, format = "BED", extraCols = extraCols)
  return(gr)
}

readCpGIslands <- function(){
  filename <- 'data/UCSCGenomeBrowser/database/cpgIslandExt.txt'
  df.cpgIslandExt <- read.table(file = filename, header = F, sep = '\t', quote = "", dec = ".", as.is = T, 
                                col.names = c('bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'length', 'cpgNum', 'gcNum', 'perCpg', 'perGc', 'obsExp'))
  gr.cpgIslandExt <- GRanges(seqnames = df.cpgIslandExt$chrom, ranges = IRanges(start = df.cpgIslandExt$chromStart, width = df.cpgIslandExt$length, names = df.cpgIslandExt$name), strand = '*',
                             bin = df.cpgIslandExt$bin, cpgNum = df.cpgIslandExt$cpgNum, gcNum = df.cpgIslandExt$gcNum, perCpg = df.cpgIslandExt$perCpg, perGc = df.cpgIslandExt$perGc, obsExp = df.cpgIslandExt$obsExp)
  return(gr.cpgIslandExt)
}

pcawg.readRNASeqData <- function(){
  x <- read.table(file = 'data/PCAWG/RNA-seq/tophat_star_fpkm_uq.v2_aliquot_gl.tsv')
}

parse.encodeDnase <- function(){
  path <- 'data/UCSCGenomeBrowser/encodeDCC/wgEncodeAwgDnaseUniform/'
  files <- list.files(path = path, pattern = '*.narrowPeak.gz')
  gr.Dnase <- foreach(i=1:length(files), .combine = c)%dopar%{
    filename <- paste0(path,files[i])
    .gr <- read.narrowPeak(filename = filename)
    return(.gr)
  }
  return(gr.Dnase)
}

kegg.getPathways <- function(organism = 'hsa'){
  pathway <- keggList(database = 'pathway', organism = organism)
  pathway_id <- substr(x = names(pathway), 6, 13)
  names(pathway) <- pathway_id
  return(pathway)
}

kegg.getGenes <- function(pathway){
  gene <- keggGet(pathway)[[1]]$GENE
  gene_id <- gene[2*(1:(length(gene)/2))-1]
  gene_name <- gene[2*(1:(length(gene)/2))]
  gene_description <- substr(x = gene_name, start = regexpr(pattern = ';', text = gene_name)+2, stop = nchar(gene_name))
  gene_name <- substr(x = gene_name, start = 1, stop = regexpr(pattern = ';', text = gene_name)-1)
  gene <- data.frame(gene_id, gene_name, gene_description)
  return(gene)
}