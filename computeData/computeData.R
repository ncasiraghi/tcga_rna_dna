library(data.table)
library(tidyverse)
library(parallel)

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

source('/BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted/computeData_functions.R')

args <- commandArgs(trailingOnly = TRUE)

sif.file <- args[1]

omic <- args[2]

dataset <- args[3]

aggregate_as  <-  args[4]

min_vaf <- as.numeric(args[5])

max_vaf <- as.numeric(args[6])

min_cov <- as.numeric(args[7])

length.sw <- as.numeric(args[8])

pos_step <- as.numeric(args[9])

custom_borders <- as.logical(args[10])

samples.in.parallel <- as.numeric(args[11])

# import and filter sif

sif <- read.delim(file = sif.file,header = TRUE,stringsAsFactors = FALSE) %>% 
  filter(data.source == dataset) %>% 
  filter(data.type == omic)

RefGenome <- unique(sif$ref.genome)

# set working dir

# outs

wd <- unique(sif$out.folder)

if(!file.exists(file.path(wd))){
  dir.create(file.path(wd), showWarnings = FALSE)
}

outdir_level_1 <- file.path(wd,paste('minvaf',min_vaf,'maxvaf',max_vaf,'mincov',min_cov,sep = '_'))

if(!file.exists(file.path(outdir_level_1))){
  dir.create(file.path(outdir_level_1), showWarnings = FALSE)
}

outdir_level_2 <- file.path(outdir_level_1,paste('aggregate',aggregate_as,sep = '_'))

if(!file.exists(file.path(outdir_level_2))){
  dir.create(file.path(outdir_level_2), showWarnings = FALSE)
}

outdir_level_3 <- file.path(outdir_level_2,paste('sw',length.sw,sep = '_'))

if(!file.exists(file.path(outdir_level_3))){
  dir.create(file.path(outdir_level_3), showWarnings = FALSE)
}

outdir_level_4 <- file.path(outdir_level_3,paste('pstep',pos_step,sep = '_'))

if(!file.exists(file.path(outdir_level_4))){
  dir.create(file.path(outdir_level_4), showWarnings = FALSE)
}

if(custom_borders){
  outdir <- file.path(outdir_level_4,paste('custom_borders'))
}else{
  outdir <- file.path(outdir_level_4,paste('default_borders'))
}

if(!file.exists(file.path(outdir))){
  dir.create(file.path(outdir), showWarnings = FALSE)
} else{
  unlink(file.path(outdir),recursive = TRUE)
  dir.create(file.path(outdir), showWarnings = FALSE)
}

setwd(outdir)

cat(paste("[",Sys.time(),"]\tOut folder:",outdir),"\n")

# pileup on snps

write.table(sif %>% select('file.snps','sample.id'),file = file.path(outdir,'samples.tsv'),sep = '\t',col.names = TRUE,row.names = FALSE,quote = FALSE)

lf <- unique(sif$file.snps)

# cytobands

if(RefGenome == 'GRCh38'){
  
  bands <- read.delim('/BCGLAB/ncasiraghi/tcga_rna_dna/data/ucsc_cytobands_GRCh38.tsv',stringsAsFactors = F,check.names = FALSE) %>% filter(chrom %in% paste0('chr',1:22))
  
}

if(RefGenome == 'GRCh37'){
  
  bands <- read.delim('/BCGLAB/ncasiraghi/tcga_rna_dna/data/ucsc_cytobands_GRCh37.tsv',stringsAsFactors = F,check.names = FALSE) %>% filter(chrom %in% paste0('chr',1:22))
  
}

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

# Define positions to aggregate data

if(custom_borders){
  
  cat(paste("[",Sys.time(),"]\tComputing CUSTOM arms borders based on available het SNPs"),"\n")
  
  listmat <- mclapply(seq_len(length(lf)),getArmsMargin,lf=lf,bands=bands,min_vaf=min_vaf,max_vaf=max_vaf,min_cov=min_cov,mc.preschedule = TRUE,mc.cores = samples.in.parallel)
  
  snps_borders <- do.call(rbind,listmat)
  
  write.table(snps_borders,file = file.path(outdir,'snps_borders.tsv'),col.names = TRUE,row.names = FALSE,sep = '\t',quote = FALSE)
  
  ab_first <- snps_borders %>%
    filter(flag == 'first') %>% 
    group_by(chr, arm) %>% 
    slice(which.max(pos))
  
  ab_last <- snps_borders %>%
    filter(flag == 'last') %>% 
    group_by(chr, arm) %>% 
    slice(which.min(pos))
  
  arm_borders <- rbind(ab_first,ab_last) %>% 
    arrange(chr,arm,pos) %>% 
    group_by(chr,arm) %>% 
    summarise(start = min(pos), end = max(pos))
  
  write.table(arm_borders,file = file.path(outdir,'arm_borders.tsv'),sep = '\t',row.names = F,quote = F,col.names = TRUE)
  
}else{
  
  cat(paste("[",Sys.time(),"]\tUsing DEFAULT arms borders based on genomic coordinates"),"\n")
  
  arm_borders <- bands %>% rename(chr = chrom)
  
}

positions <- c()

for(j in seq_len(nrow(arm_borders))){
  
  pos <- c(seq(from=arm_borders$start[j],arm_borders$end[j],by = pos_step * 1e6), arm_borders$end[j])
  
  x <- data.frame(chrom=arm_borders$chr[j],
                  arm=arm_borders$arm[j],
                  pos=unique(pos),
                  stringsAsFactors = FALSE)
  
  positions <- rbind(positions,x)
  
}

# genes of interest

# goi <- c(readLines('/BCGLAB/ncasiraghi/gene_panel/input/oncogenes.txt'),
#          readLines('/BCGLAB/ncasiraghi/gene_panel/input/cancer_genes.txt'),
#          readLines('/BCGLAB/ncasiraghi/gene_panel/input/tumor_suppressors.txt')) %>% 
#   unique() %>% 
#   sort()
# 
# ensembl <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/mart_export_GRCh38p13.tsv',check.names = F,stringsAsFactors = F,col.names = c('ensg','start','end','chr','band','Hugo_Symbol'))
# 
# ensembl <- ensembl %>% mutate(gene_center = (start + end) %/% 2)
# 
# gc <- ensembl %>% filter(Hugo_Symbol %in% goi)
# 
# gc$chr <- paste0('chr',gc$chr)
# gc$arm <- 'p'
# gc$arm[grep(gc$band,pattern = 'q')] <- 'q'
# 
# positions$dist_cancer_gene <- NA
# 
# for(idx in seq_len(nrow(positions))){
#   
#   u <- gc %>% 
#     filter(chr == positions$chrom[idx] & arm == positions$arm[idx])
#   
#   if(nrow(u) > 0){
#     positions$dist_cancer_gene[idx] <- min(abs(positions$pos[idx] - u$gene_center))
#   }
#   
# }

write.table(positions,file = file.path(outdir,'positions.tsv'),col.names = TRUE,quote = FALSE,row.names = FALSE,sep = '\t')

# run sliding window

cat(paste("[",Sys.time(),"]\tComputing data matrix:"),"\n")
cat(paste("[",Sys.time(),"]\tn. samples:\t",length(lf)),"\n")
cat(paste("[",Sys.time(),"]\tn. positions:\t",nrow(positions)),"\n")

getMatrix <- function(i,lf,bands,min_vaf,max_vaf,min_cov,length.sw,positions,aggregate_as){
  
  file <- lf[i]
  
  snps <- fread(file,data.table = FALSE,nThread = 5,verbose = FALSE) %>% 
    filter(af >= min_vaf, af <= max_vaf, cov >= min_cov)
  
  if(str_detect(snps$chr[1],'chr',negate = TRUE)){
    snps <- snps %>% 
      mutate(chr = paste0('chr',chr)) %>% 
      filter(chr %in% paste0('chr',1:22))
  }else{
    snps <- snps %>% 
      filter(chr %in% paste0('chr',1:22))
  }
  
  tomi <- which(snps$af < 0.5)
  snps$af[tomi] <- (1 - snps$af[tomi])
  
  sl <- snps %>% 
    group_by(chr) %>% 
    group_split()
  
  deck <- lapply(seq_len(length(sl)), runsw, sl, bands, length.sw = length.sw, positions = positions, aggregate_as = aggregate_as)
  
  return(as.numeric(unlist(deck)))
  
}

listmat <- mclapply(seq_len(length(lf)),getMatrix,lf=lf,bands=bands,min_vaf=min_vaf,max_vaf=max_vaf,min_cov=min_cov,length.sw=length.sw,positions=positions,aggregate_as = aggregate_as,mc.cores = samples.in.parallel)

mat <- do.call(rbind,listmat)

write.table(x = mat,file = file.path(outdir,'data.tsv'),sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)

cat(paste("[",Sys.time(),"]\t[ DONE ]"),"\n")
