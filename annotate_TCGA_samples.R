library(tidyverse)

setwd('/BCGLAB/2020_signatures/stats/annotations/')

# data from TCGA manifest
manifest <- read.delim('/BCGLAB/2020_signatures/pileup/Manifest-TCGA-BRCA.txt',stringsAsFactors = F)

manifest <- manifest %>% rename(sample = patientType, partecipant = patient)

x <- do.call(rbind,lapply(strsplit(manifest$ID,split = '-'),FUN = `[`, 1:4))

manifest$sample_vial <- apply(x,MARGIN = 1,FUN = paste,collapse='-')

manifest <- manifest %>% select(partecipant, sample, sample_vial, filename)

# downloaded files from gdc 

fl <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',pattern = '_DNA',full.names = TRUE)

file_names <- lapply(fl, function(folder) {gsub(list.files(folder,pattern = '\\.rc$'),pattern = '\\.rc',replacement = '')})

file_names <- paste0(unlist(file_names),'.bam')

manifest$bcglab <- 0
manifest$bcglab[which(manifest$filename %in% file_names)] <- 1

# missing data: in manifest but not in pileups
manifest %>% filter(bcglab == 0)

# tumor purity estimations from:
tc <- read.delim('/BCGLAB/2020_signatures/stats/annotations/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                 stringsAsFactors = F,
                 check.names = FALSE) %>% 
  select(array,purity,ploidy)

tc <- rename(tc,sample = array)

manifest <- manifest %>% left_join(.,y = tc, by = 'sample')

natcom <- read.delim('/BCGLAB/2020_signatures/stats/annotations/41467_2015_BFncomms9971_MOESM1236_ESM.txt',
                     stringsAsFactors = F,check.names = T) %>% 
  filter(Cancer.type == 'BRCA') %>% 
  rename(sample_vial = Sample.ID) %>% 
  select(-Cancer.type)

manifest <- manifest %>% left_join(.,y = natcom, by = 'sample_vial')

write.table(manifest,file = 'samples_info_file.txt',quote = F,row.names = F,col.names = T,sep = '\t')

manifest <- read.delim(file = 'samples_info_file.txt',stringsAsFactors = F)


