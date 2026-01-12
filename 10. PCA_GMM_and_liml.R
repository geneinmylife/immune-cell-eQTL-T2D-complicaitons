library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)





disease_list <- c('hypo','dka','dkd','dn','dpa','dr')

disease_info <- data.frame(upper_name=c('HYPOGLYC','KETOACIDOSIS',
                                        'NEPHROPATHY_EXMORE','NEUROPATHY',
                                        'PERIPH_ANGIOPATHY','RETINOPATHY_EXMORE'),
                           name=c('hypo','dka','dkd','dn','dpa','dr'),
                           ncase=c(7332,7841,4111,2843,2514,10413),
                           ncontrol=c(271817,271817,308539,271817,271817,308633)
)


res_data_total <- data.frame(fread('res.csv'))
IV <- read.csv('IV.csv')
R9_data_list_total <- load('Finngen_dat')





mvpcaliml_res_list <- list()
num1 <- 1
for(d1 in 1:length(disease_info$name)){
  tmp_res_d1 <- subset(res_data_total,outcome==disease_info$name[d1])
  tmp_gene_list_1 <- unique(tmp_res_d1$id)
  R9_data <- R9_data_list_total[[d1]]
  for(g1 in 1:length(tmp_gene_list_1)){
    tmp_res_d2 <- subset(tmp_res_d1,id==tmp_gene_list_1[g1])
    celltype_tmp <- unique(tmp_res_d2$celltype)
    if(length(celltype_tmp)>1){                           
      res_list_2 <- list()
      SNP_total_list <- list()
      for(c1 in 1:length(celltype_tmp)){
        tmp_gene_data <- fread('data_for_MV_IVW_PCA.csv')
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
      df_separated_tmp_summ <- subset(df_separated_tmp_summ,minP<0.05)   
      df_separated_tmp <- subset(df_separated_tmp,ID %in% df_separated_tmp_summ$SNP)
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
      if(nrow(harmonise_dat)>0){
        ld_mat <- ieugwasr::ld_matrix(unique(harmonise_dat$SNP),plink_bin = 'E:/ZJlab/software/plink/plink.exe',
                                      bfile = 'E:/ZJlab/project/DYLsproject/data/g1000_eur/g1000_eur')
      } else{
        tmp_res_d2$mvpcaliml_Pvalue <- NA
        tmp_res_d2$mvpcaliml_Beta <- NA
        tmp_res_d2$mvpcaliml_SE <- NA
        tmp_res_d2$PCGMM_Pvalue <- NA
        tmp_res_d2$PCGMM_Beta <- NA
        tmp_res_d2$PCGMM_SE <- NA
        mvpcaliml_res_list[[num1]] <- tmp_res_d2
        num1 <- num1+1
        next
      }


      beta_for_liml_PCA <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))#beta值
      se_for_liml_PCA <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      for(j in 1:length(celltype_tmp)){
        harmonise_dat_tmp <- subset(harmonise_dat,celltype==celltype_tmp[j])
        rownames(harmonise_dat_tmp) <- harmonise_dat_tmp$SNP
        harmonise_dat_tmp <- harmonise_dat_tmp[ld_mat_SNP,]
        beta_for_liml_PCA[,j] <- harmonise_dat_tmp$beta.exposure
        se_for_liml_PCA[,j] <- harmonise_dat_tmp$se.exposure
      }
      se_y <- harmonise_dat_tmp$se.outcome
      beta_y <- harmonise_dat_tmp$beta.outcome
      
      zx <- beta_for_liml_PCA/se_for_liml_PCA
      cor.x <- cor(zx)
      nx <- rep(91,ncol(beta_for_liml_PCA))
      ny <- disease_info$ncase[d1]+disease_info$ncontrol[d1]
      ld <- as.matrix(ld_mat)
      ld <- round(ld,3)
      bx <- beta_for_liml_PCA
      sx <- se_for_liml_PCA
      sy <- se_y
      
      f <- try(pca_res0 <- PCA_GMM(bx=bx,sx=sx,by=beta_y,sy=se_y,ld=ld,cor.x=cor.x,nx=nx,ny=ny,r=pca.no(0.999)),silent = T)
      if(class(f)=='try-error'){
        ci <- cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA)
        tmp_res_d2$PCGMM_Pvalue <- NA
        tmp_res_d2$PCGMM_Beta <- NA
        tmp_res_d2$PCGMM_SE <- NA
      } else{
        ci <- cbind(pca_res0$liml,pca_res0$se.liml,pca_res0$liml - qnorm(1 - 0.05/2)*pca_res0$se.liml,pca_res0$liml + qnorm(1 - 0.05/2)*pca_res0$se.liml,2*(1-pnorm(abs(pca_res0$liml/pca_res0$se.liml))),rep(0.999,ncol(bx)),rep(pca_res0$factors,ncol(bx)),rep(1-pchisq(pca_res0$Q,pca_res0$factors-(ncol(bx)+1)),ncol(bx)),pca_res0$condF,(1-pchisq((pca_res0$condF*(pca_res0$factors-ncol(bx)+1)),pca_res0$factors-ncol(bx)+1)))
        colnames(ci) <- c("estimate",'se',"lower","upper","pvalue","prop. var. expl.","PCs", "J test p-value", "cond. F-statistic","cond. F p-value")
        rownames(ci) <- tmp_res_d2$celltype
        ci <- data.frame(ci)
        tmp_res_d2$PCGMM_Pvalue <- ci$pvalue
        tmp_res_d2$PCGMM_Beta <- ci$estimate
        tmp_res_d2$PCGMM_SE <- ci$se
      }
      
      
      f1 <- try(mvpca_liml <- mv_pca_liml(Bx=bx,Sx=sx,By=beta_y,Sy=se_y,rho=ld,Phi = cor.x,R=pca.no(0.999)),silent = TRUE)
      if(class(f1)=='try-error'){
        tmp_res_d2$mvpcaliml_Pvalue <- NA
        tmp_res_d2$mvpcaliml_Beta <- NA
        tmp_res_d2$mvpcaliml_SE <- NA
        
      } else{
        mvpcaliml_res_1 <- data.frame(celltype=tmp_res_d2$celltype,
                                  Pval=mvpca_liml$pval,
                                  Beta=mvpca_liml$est,
                                  SE=mvpca_liml$se)
        tmp_res_d2$mvpcaliml_Pvalue <- mvpcaliml_res_1$Pval[match(tmp_res_d2$celltype,mvpcaliml_res_1$celltype)]
        tmp_res_d2$mvpcaliml_Beta <- mvpcaliml_res_1$Beta[match(tmp_res_d2$celltype,mvpcaliml_res_1$celltype)]
        tmp_res_d2$mvpcaliml_SE <- mvpcaliml_res_1$SE[match(tmp_res_d2$celltype,mvpcaliml_res_1$celltype)]
        
      }
      mvpcaliml_res_list[[num1]] <- tmp_res_d2
      num1 <- num1+1
    }
  }
}




























































