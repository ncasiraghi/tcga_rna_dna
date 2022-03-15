library(data.table)
library(tidyverse)

setwd('/BCGLAB/2020_signatures/synggen_inputs')

# TCGA and SPICE ids

ids <- fread('/BCGLAB/2020_signatures/synggen_inputs/resources/ids_spice_tcga_brca.txt',stringsAsFactors = F,header = F,select = c(1,3:6),na.strings = "",col.names = c('spice_id','tcga_id','clonet_purity','clonet_ploidy','kit'))

# sample info file

sif <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/input_sif_synthetic_dilutions.tsv',stringsAsFactors = F)

sif <- sif %>% mutate(plasma.snps=NA,germline.snps=NA,plasma.seg=NA,clonet.purity=NA,clonet.ploidy=NA,kit=NA,kit.bed=NA)

# seg file

seg <- fread(input = '/BCGLAB/2020_signatures/synggen_inputs/resources/brca_genomic_data.seg', 
             data.table = FALSE,stringsAsFactors = FALSE,colClasses = list(numeric=2),nThread = 4,fill = TRUE,sep = '\t') %>% distinct()

# get info

for(i in seq_len(nrow(sif))){
  
  message(sif$patient[i])
  
  info <- ids %>% filter(tcga_id == sif$patient[i])
  
  if(nrow(info) > 0){
    
    # seg file
    
    m <- seg %>% 
      filter(sample_id == unique(info$spice_id)) %>% 
      filter(segment_id != "") %>% 
      separate(segment_id,into = c('ID','chrom','loc.start','loc.end')) %>% 
      mutate(loc.start = as.numeric(loc.start),loc.end = as.numeric(loc.end)) %>% 
      select(c('ID','chrom','loc.start','loc.end','as_cn_disc','log2_corr')) %>% 
      arrange(ID,chrom,loc.start,loc.end) %>% 
      mutate(ID = sif$plasma[i])
    
    sif$plasma.seg[i] <- file.path('/BCGLAB/2020_signatures/synggen_inputs/seg',paste0(sif$plasma[i],'.seg'))
    
    write.table(m,file = sif$plasma.seg[i],col.names = TRUE,row.names = FALSE,quote = FALSE,sep = '\t')
    
    # pacbam
    
    BED <- '/BCGLAB/2020_signatures/biomart_dbsnp/mart_export_GRCh38.sort.merged.chr.bed'
    VCF <- '/BCGLAB/2020_signatures/biomart_dbsnp/dbsnp_151_all_biomart_GRCh38.chr.vcf'
    FASTA <- '/CIBIO/sharedRL/Projects/ASE/Annotations/GRCh38.d1.vd1.fa'
    OUT <- '/BCGLAB/2020_signatures/synggen_inputs/pacbam'
    PACBAM <- '/CIBIO/sharedRL/Projects/PaCBAM/git_repo/pacbam/pacbam'
    
    cmd <- paste0(PACBAM,' bam=',sif$plasma.bam[i],' bed=',BED,' vcf=',VCF,' fasta=',FASTA,' out=',OUT,' threads=20 mbq=20 mrq=20 mode=2')
    system(cmd)
    sif$plasma.snps[i] <- file.path(OUT,gsub(basename(sif$plasma.bam[i]),pattern = '\\.bam$',replacement = '.snps'))
    
    cmd <- paste0(PACBAM,' bam=',sif$germline.bam[i],' bed=',BED,' vcf=',VCF,' fasta=',FASTA,' out=',OUT,' threads=20 mbq=20 mrq=20 mode=2')
    system(cmd)
    sif$germline.snps[i] <- file.path(OUT,gsub(basename(sif$germline.bam[i]),pattern = '\\.bam$',replacement = '.snps'))
    
    # annotate 
    
    sif$clonet.ploidy[i] <- unique(info$clonet_ploidy)
    sif$clonet.purity[i] <- unique(info$clonet_purity)
    sif$kit[i] <- unique(info$kit)
    sif$kit.bed[i] <- '/BCGLAB/2020_signatures/synggen_inputs/resources/Gencode.Exome.hg38.bed'
    
  }
  
}

write.table(sif,file = 'synggen_samples_info_file.tsv',col.names = TRUE,row.names = FALSE,quote = FALSE,sep = '\t')



