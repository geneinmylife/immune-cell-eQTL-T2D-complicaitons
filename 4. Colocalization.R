library(gwasglue)
library(data.table)
library(gwasvcf)
library(biomaRt)
library(Rhtslib)
library(GenomicRanges)
library(Matrix)
library(survival)
library(calibrate)
library(ggrepel)
library(ggthemes)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(png)
library(data.table)
library(dplyr)
library(stringr)
library(rlang)
library(gwasglue)
library(usethis)
library(remotes)
library(arrow)
library(coloc)
library(tidyverse)
library("data.table")

#colocalization for dm as an example

parent_dir <- "/Users/dingyilan/Documents/result/dm"
outcome2 <- fread('/Users/dingyilan/Documents/result/dm.csv')
celltype_list <- list.dirs(parent_dir, full.names = FALSE, recursive = FALSE)
for (celltype in celltype_list) {
  setwd(paste0(parent_dir, "/", celltype))
  data_file <- paste0(celltype, ".csv") 
  if (file.exists(data_file)) {
    data_celltype <- fread(data_file, head = TRUE)
    file_list <- list.files(pattern = "*.txt")
    list <- as.data.frame(file_list)
    for (i in file_list) {
      data <- read.table(i, header = T, sep = ' ')
      colnames(data) <- c("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO")
      dat <- as.character(unlist(strsplit(data$INFO, split = ";")))
      matrix<-matrix(data=dat,ncol=6,byrow=T)
      frame<-data.frame(matrix) 
      
      dat1 <- as.character(unlist(strsplit(frame$X1, split = "=")))
      matrix1<-matrix(data=dat1,ncol=2,byrow=T)
      frame1<-data.frame(matrix1) 
      column1 <- frame1$X2
      column1<-data.frame(column1) 
      colnames(column1)[colnames(column1) == "column1"] <- "GENE"
      
      dat2 <- as.character(unlist(strsplit(frame$X2, split = "=")))
      matrix2<-matrix(data=dat2,ncol=2,byrow=T)
      frame2<-data.frame(matrix2) 
      column2 <- frame2$X2
      column2<-data.frame(column2) 
      colnames(column2)[colnames(column2) == "column2"] <- "GENESYMBOL"
      
      dat3 <- as.character(unlist(strsplit(frame$X3, split = "=")))
      matrix3<-matrix(data=dat3,ncol=2,byrow=T)
      frame3<-data.frame(matrix3) 
      column3 <- frame3$X2
      column3<-data.frame(column3) 
      colnames(column3)[colnames(column3) == "column3"] <- "PVALUE"
      
      dat4 <- as.character(unlist(strsplit(frame$X4, split = "=")))
      matrix4<-matrix(data=dat4,ncol=2,byrow=T)
      frame4<-data.frame(matrix4) 
      column4 <- frame4$X2
      column4<-data.frame(column4) 
      colnames(column4)[colnames(column4) == "column4"] <- "BETA"
      
      exp<-cbind(data,column1,column2,column3,column4)
      df <- subset(exp, select = -INFO)
      
      df$PVALUE<- as.numeric(df$PVALUE)
      df$BETA<- as.numeric(df$BETA)
      df$SE <- sqrt(((df$BETA)^2)/qchisq(df$PVALUE, 1, lower.tail = F))
      exposure <-df
      f <- gsub("\\.txt$", "", i)
      exposure1 <- subset(data_celltype, exposure == f)  
      exposure<-subset(exposure,POS>exposure1$pos.exposure-500000 & POS<exposure1$pos.exposure+500000)
      commonsnp <- outcome2[outcome2$SNP %in% exposure$ID,]
      dim(commonsnp)
      comm_snp <- as.vector(commonsnp$SNP)
      exposure <- exposure[exposure$ID %in% comm_snp,]
      outcome <- outcome2[outcome2$SNP %in% comm_snp,]
      outcome1 <- as.data.frame(outcome)
      exposure_formatted <- format_data(exposure, 
                                        type ="exposure", 
                                        snps = NULL, 
                                        header = TRUE, 
                                        snp_col = "ID",
                                        beta_col = "BETA", 
                                        se_col = "SE", 
                                        gene_col="GENESYMBOL",                                  
                                        effect_allele_col = "ALT", 
                                        other_allele_col = "REF", 
                                        pval_col = "PVALUE",                                    
                                        chr_col = "CHROM", 
                                        pos_col = "POS"
      )
      
      outcome_formatted <- format_data(outcome1,type ="outcome",  snps = NULL, header = TRUE, snp_col = "SNP",beta_col = "Beta", se_col = "Se", effect_allele_col = "Effect_allele",other_allele_col = "Other_allele",pval_col = "P.x",eaf_col="Effect_allele_freq" )
      dat<-harmonise_data(exposure_dat = exposure_formatted,outcome_dat = outcome_formatted,action=1)
      dat$MAF <- ifelse (dat$eaf.outcome<=0.5, as.numeric(dat$eaf.outcome),
                         ifelse(dat$eaf.outcome>=0.5, as.numeric(1-(dat$eaf.outcome)), NA))
      coloc_results <- coloc.abf(dataset2=list(snp=dat$SNP,pvalues=dat$pval.outcome, type="cc", s=0.133709, N=1812017,MAF=dat$MAF),
                                 dataset1=list(snp=dat$SNP,pvalues=dat$pval.exposure, type="quant", N=91,MAF=dat$MAF))
      coloc_results$summary
      head(coloc_results$results)
      coloc_results_summary <- as.data.frame(coloc_results$summary)
      coloc_results_summary
      coloc_results_all <- as.data.frame(coloc_results$results)
      
      dir <- file.path( "result")
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        cat("Created directory:", dir, "\n")
      } 
      
      
      dir_path <- paste0(dir, "/")
      
      if (!dir.exists(dir_path)) {
        dir.create(dir_path)
      }
      
      filename1 <- file.path(dir, paste("coloc_results_", f, "_final.txt", sep = "")) 
      filename2 <- file.path(dir, paste("coloc_results_all_", f, ".txt", sep = "")) 
      write.table(coloc_results_summary, filename1, quote=F, row.names=T, sep="\t")
      write.table(coloc_results_all, filename2, quote=F, row.names=F, sep="\t")
      
      setwd(paste0(parent_dir, "/", celltype))
    }
    
  } else {
    cat("File not found:", data_file, "\n")
  }
}
