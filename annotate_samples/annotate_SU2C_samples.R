library(dplyr)

setwd('/BCGLAB/2020_signatures/annotations/SU2C/')

m <- read.delim('/BCGLAB/2020_signatures/annotations/SU2C/150_cases_SU2C_PMID_26000489_TabS2.txt') %>% rename(patient.id = PATIENT_ID)

sif.file <- '/BCGLAB/ncasiraghi/tcga_rna_dna/data/datasetSNPs.tsv'

sif <- read.delim(file = sif.file,header = TRUE,stringsAsFactors = FALSE)

# SU2C

ctrls <- readLines('/BCGLAB/2020_signatures/pileup/SU2C/BAMs_Normal.txt')
cases <- readLines('/BCGLAB/2020_signatures/pileup/SU2C/BAMs_Tumor.txt')

su2c <- sif %>% filter(data.source == 'SU2C')

su2c$is.tumor <- 0
su2c$is.tumor[which(su2c$sample.id %in% gsub(basename(cases),pattern = '\\.bam',replacement = ''))] <- 1

su2c$is.tissue <- 1

su2c <- su2c %>% select(sample.id,is.tumor,is.tissue,data.source)

su2c$patient.id <- NA

idx <- grep(su2c$sample.id,pattern = '_')
su2c$patient.id[idx] <- sapply(lapply(strsplit(su2c$sample.id[idx],split = '_'),FUN = `[`,1:2),FUN = paste,collapse = '_')

idx <- grep(su2c$sample.id,pattern = '_',invert = TRUE)
su2c$patient.id[idx] <- as.numeric(sapply(strsplit(su2c$sample.id[idx],split = '-'),FUN = `[`,1))

su2c <- left_join(x = su2c,y = m,by='patient.id')

write.table(su2c,file = 'SU2C_sif.tsv',quote = FALSE,sep = '\t',row.names = F,col.names = TRUE)