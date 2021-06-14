library(tidyverse)

setwd('/BCGLAB/2020_signatures/stats/annotations/TCGA-BRCA/')

# data from TCGA manifest
original_manifest <- read.delim('/BCGLAB/2020_signatures/manifest.txt',stringsAsFactors = F)

manifest <- original_manifest %>% 
  rename(TCGA_barcode = ID.x, sample = patientType, partecipant = patient, filename = file) %>% 
  select(-ID.y, -Patient)

x <- do.call(rbind,lapply(strsplit(manifest$TCGA_barcode,split = '-'),FUN = `[`, 1:4))

manifest$sample_vial <- apply(x,MARGIN = 1,FUN = paste,collapse='-')

manifest <- manifest %>% filter(Disease == 'TCGA-BRCA')

# add subtype from cbioportal

tt <- read.delim('/BCGLAB/mimesis/panel/cbioportal/outs/sif_cbioportal_brca.tsv',stringsAsFactors = F) %>% 
  select(data,sample.id,class,type) %>% 
  rename(sample = sample.id, tumor_type = type)

manifest <- left_join(manifest,tt,by = 'sample')
  
# downloaded files from gdc 

fl <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',pattern = '_DNA',full.names = TRUE)

file_names <- lapply(fl, function(folder) {gsub(list.files(folder,pattern = '\\.rc$'),pattern = '\\.rc',replacement = '')})

file_names <- paste0(unlist(file_names),'.bam')

manifest$bcglab <- 0
manifest$bcglab[which(manifest$filename %in% file_names)] <- 1

miss <- manifest %>% filter(bcglab == 0) %>% pull(filename)

write.table(original_manifest %>% filter(file %in% miss), file = 'BRCA_missing_samples.tsv',quote = F,row.names = F,col.names = T,sep = '\t')

# tumor purity estimations from: http://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5
tc <- read.delim('/BCGLAB/2020_signatures/stats/annotations/TCGA-BRCA/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                 stringsAsFactors = F,
                 check.names = FALSE) %>% 
  select(array,purity,ploidy)

tc <- rename(tc,sample = array)

manifest <- manifest %>% left_join(.,y = tc, by = 'sample')

# tumor purity estimations from: https://www.nature.com/articles/ncomms9971
natcom <- read.delim('/BCGLAB/2020_signatures/stats/annotations/TCGA-BRCA/41467_2015_BFncomms9971_MOESM1236_ESM.txt',
                     stringsAsFactors = F,check.names = T) %>% 
  filter(Cancer.type == 'BRCA') %>% 
  rename(sample_vial = Sample.ID) %>% 
  select(-Cancer.type)

manifest <- manifest %>% 
  left_join(.,y = natcom, by = 'sample_vial') %>% 
  rename(TCGA_purity = purity, TCGA_ploidy = ploidy)
  
# studio tumor purity and selection

# batch effect affected samples

to.remove <- c('TCGA-3C-AAAU-01_DNA_tumor', 'TCGA-3C-AAAU-10_DNA_normal',
               'TCGA-3C-AALI-01_DNA_tumor', 'TCGA-3C-AALI-10_DNA_normal',
               'TCGA-3C-AALJ-01_DNA_tumor', 'TCGA-3C-AALJ-10_DNA_normal',
               'TCGA-3C-AALK-01_DNA_tumor', 'TCGA-3C-AALK-10_DNA_normal',
               'TCGA-4H-AAAK-01_DNA_tumor', 'TCGA-4H-AAAK-10_DNA_normal',
               'TCGA-5L-AAT0-01_DNA_tumor', 'TCGA-5L-AAT0-10_DNA_normal',
               'TCGA-5L-AAT1-01_DNA_tumor', 'TCGA-5L-AAT1-10_DNA_normal',
               'TCGA-5T-A9QA-01_DNA_tumor', 'TCGA-5T-A9QA-10_DNA_normal',
               'TCGA-GI-A2C8-01_DNA_tumor', 'TCGA-GI-A2C8-11_DNA_normal',
               'TCGA-GI-A2C9-01_DNA_tumor', 'TCGA-GI-A2C9-11_DNA_normal',
               'TCGA-HN-A2NL-01_DNA_tumor', 'TCGA-HN-A2NL-10_DNA_normal',
               'TCGA-HN-A2OB-01_DNA_tumor', 'TCGA-HN-A2OB-10_DNA_normal',
               'TCGA-JL-A3YW-01_DNA_tumor', 'TCGA-JL-A3YW-10_DNA_normal',
               'TCGA-JL-A3YX-01_DNA_tumor', 'TCGA-JL-A3YX-10_DNA_normal',
               'TCGA-LQ-A4E4-01_DNA_tumor', 'TCGA-LQ-A4E4-10_DNA_normal',
               'TCGA-MS-A51U-01_DNA_tumor', 'TCGA-MS-A51U-10_DNA_normal',
               'TCGA-OK-A5Q2-01_DNA_tumor', 'TCGA-OK-A5Q2-10_DNA_normal',
               'TCGA-PE-A5DC-01_DNA_tumor', 'TCGA-PE-A5DC-10_DNA_normal',
               'TCGA-PE-A5DD-01_DNA_tumor', 'TCGA-PE-A5DD-10_DNA_normal',
               'TCGA-PE-A5DE-01_DNA_tumor', 'TCGA-PE-A5DE-10_DNA_normal',
               'TCGA-PL-A8LV-01_DNA_tumor', 'TCGA-PL-A8LV-10_DNA_normal',
               'TCGA-PL-A8LX-01_DNA_tumor', 'TCGA-PL-A8LX-10_DNA_normal',
               'TCGA-PL-A8LY-01_DNA_tumor', 'TCGA-PL-A8LY-10_DNA_normal',
               'TCGA-PL-A8LZ-01_DNA_tumor', 'TCGA-PL-A8LZ-10_DNA_normal',
               'TCGA-UL-AAZ6-01_DNA_tumor', 'TCGA-UL-AAZ6-10_DNA_normal',
               'TCGA-UU-A93S-01_DNA_tumor', 'TCGA-UU-A93S-10_DNA_normal',
               'TCGA-V7-A7HQ-01_DNA_tumor', 'TCGA-V7-A7HQ-10_DNA_normal',
               'TCGA-W8-A86G-01_DNA_tumor', 'TCGA-W8-A86G-10_DNA_normal',
               'TCGA-WT-AB41-01_DNA_tumor', 'TCGA-WT-AB41-10_DNA_normal',
               'TCGA-WT-AB44-01_DNA_tumor', 'TCGA-WT-AB44-10_DNA_normal',
               'TCGA-XX-A899-01_DNA_tumor', 'TCGA-XX-A899-10_DNA_normal',
               'TCGA-XX-A899-11_DNA_normal', 'TCGA-XX-A89A-01_DNA_tumor',
               'TCGA-XX-A89A-10_DNA_normal', 'TCGA-Z7-A8R5-01_DNA_tumor',
               'TCGA-Z7-A8R5-10_DNA_normal', 'TCGA-Z7-A8R6-01_DNA_tumor',
               'TCGA-Z7-A8R6-10_DNA_normal')

to.remove <- gsub(to.remove,pattern = '_DNA_tumor|_DNA_normal',replacement = '')

manifest$batch_approved <- 1
manifest$batch_approved[which(manifest$sample %in% to.remove)] <- 0

manifest %>%
  filter(bcglab == 1 & batch_approved == 1) %>%
  group_by(type) %>% 
  summarise(n = n())

cnt <- manifest %>%
  filter(bcglab == 1 & batch_approved == 1 & type == 1) %>% 
  gather(method,tumor_purity,c(TCGA_purity,ESTIMATE:CPE))

p <- ggplot(cnt, aes(x=method, y=tumor_purity)) + 
  geom_violin() +
  geom_boxplot(width=0.1)

ggsave(filename = 'BRCA_tumor_purity.pdf',plot = p, width = 300,height = 150,dpi = 300,units = 'mm',device = 'pdf')

cnt %>%
  group_by(method) %>% 
  summarise(n = sum(!is.na(tumor_purity))) %>% 
  arrange(desc(n))

selection <- cnt %>%
  group_by(filename) %>% 
  summarise(n = sum(!is.na(tumor_purity)),
            mean_tumor_purity = mean(tumor_purity,na.rm = TRUE),
            sd_tumor_purity = sd(tumor_purity,na.rm = TRUE)) %>% 
  filter(n == 6) %>% 
  filter(mean_tumor_purity > 0.85, sd_tumor_purity < 0.05)

sel <- manifest %>%
  filter(bcglab == 1 & batch_approved == 1 & type == 1) %>% 
  filter(filename %in% selection$filename)

not_sel <- manifest %>%
  filter(bcglab == 1 & batch_approved == 1 & type == 1) %>% 
  filter(!filename %in% selection$filename)

compute_fisher <- function(sel,not_sel,brca_subtype){
  
  colA <- c(nrow(sel %>% filter(tumor_type == brca_subtype)), nrow(not_sel %>% filter(tumor_type == brca_subtype)))
  colB <- c(nrow(sel %>% filter(!is.na(tumor_type) & type != brca_subtype)), nrow(not_sel %>% filter(!is.na(tumor_type) & type != brca_subtype)))
  
  dat <- data.frame(brca_subtype_yes = colA,
                    brca_subtype_not = colB,
                    row.names = c('TC selected','TC not selected'),
                    stringsAsFactors = FALSE,check.names = FALSE)
  
  return(list(brca_subtype=brca_subtype,dat=dat,ftest=fisher.test(dat)))
}

compute_fisher(sel,not_sel,brca_subtype = 'HER2+')
compute_fisher(sel,not_sel,brca_subtype = 'HR+')
compute_fisher(sel,not_sel,brca_subtype = 'TNBC')

download.tumor <- original_manifest %>% filter(file %in% selection$filename)

download.match <- original_manifest %>% filter(Patient %in% unique(download.tumor$Patient) & type %in% c(10,11))
                  
write.table(rbind(download.tumor,download.match) %>% arrange(Patient),file = 'BRCA_samples_tumor_purity_selected.txt',quote = F,row.names = F,col.names = T,sep = '\t')

# update manifest file
manifest$note <- NA

# selected for dilution

manifest$note[which(manifest$filename %in% sel$filename)] <- 'dilution'
keep <- manifest %>% filter(note == 'dilution') %>% pull(partecipant)
manifest$note[which(manifest$partecipant %in% keep & manifest$type %in% c(10,11))] <- 'dilution' 

# missing tumor purity
manifest$note[which(is.na(manifest$TCGA_purity))] <- 'missing_tumor_purity'
keep <- manifest %>% filter(note == 'missing_tumor_purity') %>% pull(partecipant)
manifest$note[which(manifest$partecipant %in% keep & manifest$type %in% c(10,11))] <- 'missing_tumor_purity' 

# low tumor purity
manifest$note[which(manifest$TCGA_purity < 0.25)] <- 'low_tumor_purity'
keep <- manifest %>% filter(note == 'low_tumor_purity') %>% pull(partecipant)
manifest$note[which(manifest$partecipant %in% keep & manifest$type %in% c(10,11))] <- 'low_tumor_purity' 

manifest$sample_bcglab <- paste0(manifest$sample,'_DNA_tumor')
manifest <- relocate(manifest,sample_bcglab,.after = sample)

write.table(manifest,file = 'BRCA_samples_info_file.txt',quote = F,row.names = F,col.names = T,sep = '\t')
