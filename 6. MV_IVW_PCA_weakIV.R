library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)



#MV-IVW-PCA

# Exposure and outcome info
DICE_data_list <- list.files('DICE_unfiltered/')
celltype_list <- gsub('.vcf.gz','',DICE_data_list)
finngen_data_list <- list.files('finngen_R9_DM/data/') 
disease_list <- c('hypo','dka','dkd','dn','dpa','dr')

disease_info <- data.frame(upper_name=c('HYPOGLYC','KETOACIDOSIS',
                                        'NEPHROPATHY_EXMORE','NEUROPATHY',
                                        'PERIPH_ANGIOPATHY','RETINOPATHY_EXMORE'),
                           name=c('hypo','dka','dkd','dn','dpa','dr'),
                           ncase=c(7332,7841,4111,2843,2514,10413),
                           ncontrol=c(271817,271817,308539,271817,271817,308633)
)

#read main res 
res_data_total <- read.csv('res_data_total.csv')





mvpca_res_list <- list()
num1 <- 1
for(d1 in 1:length(disease_info$name)){
  tmp_res_d1 <- subset(res_data_total,outcome==disease_info$name[d1])
  tmp_gene_list_1 <- unique(tmp_res_d1$id)
  R9_data <- R9_data_list_total[[d1]]
  for(g1 in 1:length(tmp_gene_list_1)){
    tmp_res_d2 <- subset(tmp_res_d1,id==tmp_gene_list_1[g1])
    IV_tmp <- subset(IV,id == tmp_gene_list_1[g1])
    celltype_tmp <- unique(IV_tmp$celltype)
    if(length(celltype_tmp)>1){                           #only pleiotropy
      res_list_2 <- list()
      SNP_total_list <- list()
      for(c1 in 1:length(celltype_tmp)){
        tmp_gene_data <- fread(paste0('MV_IVW_PCA/SNP_celltype/',celltype_tmp[c1],'_data_for_MV_IVW_PCA.csv'))
        df_separated <- tmp_gene_data
        SNP_total_list[[c1]] <- df_separated
      }
      SNP_total_dat <- do.call(rbind,SNP_total_list)
      df_separated_tmp <- subset(SNP_total_dat,Gene==tmp_gene_list_1[g1])
      tmp1 <- data.frame(table(df_separated_tmp$ID))
      tmp2 <- subset(tmp1,Freq>=length(celltype_tmp))
      df_separated_tmp <- subset(df_separated_tmp,ID %in% tmp2$Var1)
      df_separated_tmp_summ <- df_separated_tmp %>% group_by(ID) %>% summarise(min(Pvalue))
      colnames(df_separated_tmp_summ) <- c('SNP','minP')
      df_separated_tmp <- subset(df_separated_tmp,ID %in% df_separated_tmp_summ$SNP)
      
      
      #Keep SNPs present in both FinnGen and DICE
      df_separated_tmp$unique_SNPID <- paste(df_separated_tmp$ID,df_separated_tmp$REF,df_separated_tmp$ALT,sep = '_')
      df_separated_tmp$unique_SNPID2 <- paste(df_separated_tmp$ID,df_separated_tmp$ALT,df_separated_tmp$REF,sep = '_')
      df_separated_tmp <- subset(df_separated_tmp,unique_SNPID %in% R9_data$unique_SNPID | unique_SNPID %in% R9_data$unique_SNPID2)
      df_separated_tmp <- subset(df_separated_tmp,ID %in% df_separated_tmp_summ$SNP)
      SNP_ID_tmp1 <- unique(df_separated_tmp$ID)
      
      df_separated_tmp <- data.frame(df_separated_tmp)
      exposure_formatted <- format_data(df_separated_tmp,
                                        type ="exposure",
                                        snps = NULL,
                                        header = TRUE,
                                        phenotype_col = "phenotype",
                                        snp_col = "ID",
                                        beta_col = "Beta",
                                        se_col = "SE",
                                        gene_col="GeneSymbol",
                                        effect_allele_col = "ALT",
                                        other_allele_col = "REF",
                                        pval_col = "Pvalue",
                                        chr_col = "CHROM",
                                        pos_col = "POS"
      )
      
      R9_dat_tmp_1 <- subset(R9_data,rsids %in% exposure_formatted$SNP)
      outcome_formatted <- format_data(R9_dat_tmp_1,
                                       type ="outcome",
                                       # phenotype='Phenotype',
                                       snps = NULL,
                                       header = TRUE,
                                       snp_col = "rsids",
                                       beta_col = "beta",
                                       se_col = "sebeta",
                                       effect_allele_col = "alt",
                                       other_allele_col = "ref",
                                       pval_col = "pval",
                                       eaf_col="af_alt",
                                       samplesize_col = 'n',
                                       phenotype_col = 'disease'
      )
      harmonise_dat<-harmonise_data(exposure_dat = exposure_formatted,outcome_dat = outcome_formatted)
      harmonise_dat <- subset(harmonise_dat,mr_keep==TRUE)
      harmonise_dat_tmpx1 <- data.frame(table(harmonise_dat$SNP))
      harmonise_dat_tmpx1 <- subset(harmonise_dat_tmpx1,Freq==length(celltype_tmp))
      harmonise_dat <- subset(harmonise_dat,SNP %in% harmonise_dat_tmpx1$Var1)
      #calculate ld matrix
      ld_mat <- ieugwasr::ld_matrix(unique(harmonise_dat$SNP),plink_bin = 'software/plink/plink.exe',
                                    bfile = 'data/g1000_eur/g1000_eur')
      
      
      
      ld_mat_SNP <- str_split_fixed(rownames(ld_mat),'_',2)[,1]
      harmonise_dat <- subset(harmonise_dat,SNP %in% ld_mat_SNP)
      harmonise_dat$celltype <- str_split_fixed(harmonise_dat$exposure,'_',2)[,2]
      

      
      
      #MV_IVW_PCA
      beta_for_MVMR_PCA <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      se_for_MVMR_PCA <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      for(j in 1:length(celltype_tmp)){
        harmonise_dat_tmp <- subset(harmonise_dat,celltype==celltype_tmp[j])
        rownames(harmonise_dat_tmp) <- harmonise_dat_tmp$SNP
        harmonise_dat_tmp <- harmonise_dat_tmp[ld_mat_SNP,]
        beta_for_MVMR_PCA[,j] <- harmonise_dat_tmp$beta.exposure
        se_for_MVMR_PCA[,j] <- harmonise_dat_tmp$se.exposure
      }
      se_y <- harmonise_dat_tmp$se.outcome
      beta_y <- harmonise_dat_tmp$beta.outcome
      beta_for_MVMR_PCA_2 <- abs(beta_for_MVMR_PCA)
      Psi = (rowSums(beta_for_MVMR_PCA_2)/se_y)%o%
        (rowSums(beta_for_MVMR_PCA_2)/se_y)*ld_mat
      
      if(length(celltype_tmp) > nrow(ld_mat)){
        tmp_res_d2$mvpca_Pvalue <- NA
        tmp_res_d2$mvpca_Beta <- NA
        tmp_res_d2$mvpca_SE <- NA
        mvpca_res_list[[num1]] <- tmp_res_d2
        num1 <- num1+1
        next
      }
      
    
      
      K2 =which(cumsum(prcomp(Psi, scale=FALSE)$sdev^2/sum((prcomp(Psi, scale=FALSE)$sdev^2)))>0.999)[1]
      
      
      beta_XG_for_MVMR_PCA <- data.frame(matrix(0,nrow = K2,ncol = length(celltype_tmp)))
      for(l in 1:length(celltype_tmp)){#
        beta_XG_for_MVMR_PCA[,l] <- as.numeric(beta_for_MVMR_PCA[,l]%*%prcomp(Psi, scale=FALSE)$rotation[,1:K2])
      }
      betaYG0 = as.numeric(beta_y%*%prcomp(Psi, scale=FALSE)$rotation[,1:K2])
      
      sebetaYG0 = as.numeric(se_y%*%prcomp(Psi, scale=FALSE)$rotation[,1:K2])
      Sigma = se_y%o%se_y*ld_mat
      pcSigma = t(prcomp(Psi, scale=FALSE)$rotation[,1:K2])%*%Sigma%*%
        prcomp(Psi, scale=FALSE)$rotation[,1:K2]
      
      f <- try(mvpca <- mr_mvivw(mr_mvinput(as.matrix(beta_XG_for_MVMR_PCA),
                                            matrix(1,nrow = K2,ncol = length(celltype_tmp)),
                                            betaYG0, rep(1, nrow(beta_XG_for_MVMR_PCA)), corr=pcSigma), model="fixed"),silent = TRUE)
      if(class(f)=='try-error'){
        tmp_res_d2$mvpca_Pvalue <- NA
        tmp_res_d2$mvpca_Beta <- NA
        tmp_res_d2$mvpca_SE <- NA
        
      } else{
        mvpca_res_1 <- data.frame(celltype=celltype_tmp,
                                  Pval=mvpca$Pvalue,
                                  Beta=mvpca$Estimate,
                                  SE=mvpca$StdError)
        tmp_res_d2$mvpca_Pvalue <- mvpca_res_1$Pval[match(tmp_res_d2$celltype,mvpca_res_1$celltype)]
        tmp_res_d2$mvpca_Beta <- mvpca_res_1$Beta[match(tmp_res_d2$celltype,mvpca_res_1$celltype)]
        tmp_res_d2$mvpca_SE <- mvpca_res_1$SE[match(tmp_res_d2$celltype,mvpca_res_1$celltype)]
        
      }
      
      
      mvpca_res_list[[num1]] <- tmp_res_d2
      num1 <- num1+1
    }
  }
}
