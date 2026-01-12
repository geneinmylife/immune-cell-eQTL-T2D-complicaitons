library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(ieugwasr)
library(data.table)

setwd("/Users/dingyilan/Documents/DICE")
dice <- fread("dice.csv")
dice_clean <- dice %>% filter((chr.exposure != 'chr6' ) | (chr.exposure== 'chr6' & (pos.exposure < 26000000 |pos.exposure > 34000000))) 
write.csv(dice_clean,"dice_clean.csv")
exposure_formatted_clump <- fread("dice_clean.csv")
#dm
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/T2DM.csv",sep = ",",gene_col = "gene.exposure",snp_col = "SNP",beta_col = "Beta",se_col = "Se",eaf_col="Effect_allele_freq",effect_allele_col = "Effect_allele",other_allele_col = "Other_allele",pval_col = "P.x")
#dpa
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/dpa/finngen_R9_DM_PERIPH_ANGIOPATHY.csv",sep = ",",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",eaf_col="af_alt",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval")
#dn
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/dn/finngen_R9_DM_NEUROPATHY.csv",sep = ",",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",eaf_col="af_alt",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval")
#dkd
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/dkd/finngen_R9_DM_NEPHROPATHY_EXMORE.csv",sep = ",",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",eaf_col="af_alt",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval")
#dr
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/dr/finngen_R9_DM_RETINOPATHY_EXMORE.csv",sep = ",",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",eaf_col="af_alt",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval")
#dka
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/dka/finngen_R9_DM_KETOACIDOSIS.csv",sep = ",",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",eaf_col="af_alt",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval")
#hypo
outcome_dat<-read_outcome_data(snps = exposure_formatted_clump$SNP,filename = "/Users/dingyilan/Documents/hypo/finngen_R9_DM_HYPOGLYC.csv",sep = ",",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",eaf_col="af_alt",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval")

dat<-harmonise_data(exposure_dat = exposure_formatted_clump,outcome_dat = outcome_dat)
res <- mr(dat)
res
OR <-generate_odds_ratios(res)
res_single <- mr_singlesnp(dat)
ORR <-generate_odds_ratios(res_single)

setwd('/Users/dingyilan/Documents/result')
file_list <- list.files(pattern = "*.csv")
list <- as.data.frame(file_list)
result <- NULL
for (i in file_list) {
  data <- read.csv(i)
  data$outcome <- i
  result <- rbind(result,data)
} 
result$merge <- paste(result$exposure,result$method,result$outcome,sep = "_")
data <- result[,c(8:19,21,22)]
data <- data[!duplicated(data), ]
merge <- result[,c(2:7,22)]
data$fdr_new <- p.adjust(data$pval,method = 'BH')
dat <- merge(merge,data,by.x='merge',by.y='merge')
dat_fdr <- subset(dat,fdr_new<0.05)
dat_fdr$merge <- paste(dat_fdr$exposure,dat_fdr$method,sep = "_")
dat_fdr <- dat_fdr[,c(1:21)]
write.csv(dat,'result_fdr.csv')
