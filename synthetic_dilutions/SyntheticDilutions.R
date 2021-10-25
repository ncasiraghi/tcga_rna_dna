library(parallel)
library(data.table)
library(tidyverse)

# TCGA-BRCA selected samples

setwd('/BCGLAB/2020_signatures/synthetic_dilutions')

brca <- read.delim('/BCGLAB/2020_signatures/annotations/TCGA-BRCA/BRCA_samples_info_file.txt',stringsAsFactors = FALSE)

brca_tumor <- brca %>% 
  filter(note == 'dilution', batch_approved == 1, type == 1) %>% 
  select(partecipant,filename,sample_bcglab,TCGA_purity,id)

brca_normal <- brca %>% 
  filter(note == 'dilution', batch_approved == 1, type != 1) %>% 
  select(partecipant,filename,sample_bcglab,id)

m <- full_join(x = brca_tumor,y = brca_normal, by = 'partecipant',suffix = c('_case','_ctrl'))

bamfolder <- '/BCGLAB/2020_signatures/synthetic_dilutions/bams'

df <- data.frame(patient = m$partecipant,
           plasma = m$sample_bcglab_case,
           plasma.bam = file.path(bamfolder,m$id_case,m$filename_case),
           germline = m$sample_bcglab_ctrl,
           germline.bam = file.path(bamfolder,m$id_ctrl,m$filename_ctrl),
           germline.mean.coverage = NA,
           plasma.mean.coverage = NA,
           TC = m$TCGA_purity,
           stringsAsFactors = FALSE)

head(df)

# load coverage

cov <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/coverage/snps_coverage_per_sample.tsv',stringsAsFactors = FALSE)

for(i in seq_len(nrow(df))){
  
  cov.case <- cov$mean_cov[which(cov$sample == df$plasma[i])]
  cov.germ <- cov$mean_cov[which(cov$sample == df$germline[i])]
  
  df$plasma.mean.coverage[i] <- cov.case
  df$germline.mean.coverage[i] <- cov.germ
  
}

sif <- df %>% 
  arrange(desc(germline.mean.coverage),plasma.mean.coverage,TC)

write.table(sif,file = 'input_sif_synthetic_dilutions.tsv',col.names = TRUE,row.names = FALSE,quote = FALSE,sep = '\t')

compute_adm <- function(tc){
  adm <- as.character(1-tc)
  if(nchar(adm) == 3){
    adm <- paste0(substr(adm, 3, 3), "0", "_0")
  } else if (nchar(adm) == 4){
    adm <- paste0(substr(adm, 3, 4), "_0")
  } else if (nchar(adm) == 5){
    adm <- paste0(substr(adm, 3, 4), "_", substr(adm, 5, 5))
  }
  return(adm)
} #takes in input a tumor content (0<tc<1) and returns the admixture in "xx_x" format

subsample_bam <- function(i, adm){
  if(!is.na(p.plasma[i]) && !is.na(p.germ[i]) && p.plasma[i] != 1 && p.germ[i] < 1){
    message(paste0(Sys.time(), "\tComputing Admixture: ", adm))
    samtools_step1_plasma <- paste0("samtools view -s ", p.plasma[i], " -b ", sif$plasma.bam[i], " > ",
                                    file.path(path_bam, patient_folder[i]), "/to_merge_", plasma_bam_name[i])
    samtools_step1_germline <- paste0("samtools view -s ", p.germ[i], " -b ", sif$germline.bam[i], " > ",
                                      file.path(path_bam, patient_folder[i]), "/to_merge_GERM_", plasma_bam_name[i])
    samtools_step2 <- paste0("samtools merge ", file.path(path_bam, patient_folder[i]), "/", 
                             substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".bam ",
                             file.path(path_bam, patient_folder[i]), "/to_merge_", plasma_bam_name[i], " ",
                             file.path(path_bam, patient_folder[i]), "/to_merge_GERM_", plasma_bam_name[i])
    samtools_sort <- paste0("samtools sort ", file.path(path_bam, patient_folder[i]), "/", 
                            substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".bam ", #samtools v1.19 requires >
                            "-o ",file.path(path_bam, patient_folder[i]), "/", 
                            substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".sorted.bam") #samtools v1.19 requires .bam
    rm_bam1 <- paste0("rm ", file.path(path_bam, patient_folder[i]), "/to_merge_", plasma_bam_name[i])
    rm_bam2 <- paste0("rm ", file.path(path_bam, patient_folder[i]), "/to_merge_GERM_", plasma_bam_name[i])
    rm_bam3 <- paste0("rm ", file.path(path_bam, patient_folder[i]), "/", 
                      substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".bam ")
    samtools_index <- paste0("samtools index ", file.path(path_bam, patient_folder[i]), "/", 
                             substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".sorted.bam")
    
    
    if(FALSE){
      pacbam_adm <- paste0("/elaborazioni/sharedCO/CO_Shares/Code/PaCBAM/pacbam bam=", path_bam, patient_folder[i], "/",
                           substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".sorted.bam ",
                           "bed=/CIBIO/sharedCO/PaCBAM/PCF_SELECT/Nimblegen_Regions.bed vcf=/CIBIO/sharedCO/PaCBAM/PCF_SELECT/Nimblegen_Regions_nodup.vcf fasta=/elaborazioni/sharedCO/PCF_Project/Alignment_Pipeline/human_g1k_v37.fasta strandbias mode=0 threads=5 mbq=20 mrq=20 regionperc=0.5 out=", path_bam, "PaCBAM/")
      
      # CREATE SIF FILE FOR NEXT STEPS
      write.table(cbind(paste0(path_bam, "PaCBAM/", substr(plasma_bam_name[i], 1, nchar(plasma_bam_name[i])-4), "_adm", adm, ".sorted.snps"),
                        paste0("/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/PaCBAM/Cornell/", gsub(".bam", ".snps", germline_bam_name[i]))), 
                  file = "/SPICE/downloads/PCF_SELECT/SyntheticDilutions/SIF.tsv", append = T, sep = "\t", quote = F, col.names = F, row.names = F)
    }
    
    system(samtools_step1_plasma)
    system(samtools_step1_germline)
    system(samtools_step2)
    system(samtools_sort)
    system(rm_bam1)
    system(rm_bam2)
    system(rm_bam3)
    system(samtools_index)
    # system(pacbam_adm)
    # system(pacbam)
  }
}
#extract strings for samtools
path_bam <- "/BCGLAB/2020_signatures/synthetic_dilutions/insilico_data"
# patient_folder <- do.call(c, lapply(1:nrow(sif), function(x) strsplit(sif$plasma.bam, "/")[[x]][8]))
patient_folder <- sif$patient
plasma_bam_name <- basename(sif$plasma.bam)
germline_bam_name <- basename(sif$germline.bam)

if(TRUE){
  for(name in patient_folder){
    dir.create(file.path(path_bam, name),showWarnings = FALSE) 
  }
}

cores = 20

tc.in <- sif$TC
# TC <- seq(0.01,0.2,0.02)
TC <- c(0.10,0.05,0.01)

for(tc in TC){
  # Compute subsampling percentage for plasma and germline to obtain 
  # a new bam with the desired Tumor Content
  tc.fin <- rep(tc, nrow(sif))
  germline.mean.cov <- sif$germline.mean.coverage
  plasma.mean.cov <- sif$plasma.mean.coverage
  # plasma.mean.cov : tc.in*plasma.mean.cov = X : tc.fin*plasma.mean.cov
  plasma.reads.fin <- (plasma.mean.cov*(tc.fin*plasma.mean.cov))/(tc.in*plasma.mean.cov)
  perc.samtools.plasma <- plasma.reads.fin/plasma.mean.cov
  perc.samtools.germ <- ((1-perc.samtools.plasma)*plasma.mean.cov)/germline.mean.cov
  
  p.plasma <- ifelse(perc.samtools.plasma > 1 | perc.samtools.plasma < 0, NA, perc.samtools.plasma)
  p.germ <- ifelse(perc.samtools.germ > 1 | perc.samtools.germ < 0, NA, perc.samtools.germ)
  ### perc.samtools.germ > 1 --- remove
  
  adm <- compute_adm(tc)
  
  mclapply(seq_len(length(p.plasma)), subsample_bam, adm = adm, mc.cores = cores)
}
