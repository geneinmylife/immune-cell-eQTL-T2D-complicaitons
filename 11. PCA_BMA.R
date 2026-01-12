library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)


mvpcaliml_res_list <- list()
num1 <- 1
for(d1 in 1:length(disease_info$name)){#
  tmp_res_d1 <- subset(res_data_total,outcome==disease_info$name[d1])
  tmp_gene_list_1 <- unique(tmp_res_d1$id)
  R9_data <- R9_data_list_total[[d1]]
  for(g1 in 1:length(tmp_gene_list_1)){#
    tmp_res_d2 <- subset(tmp_res_d1,id==tmp_gene_list_1[g1])
    celltype_tmp <- unique(tmp_res_d2$celltype)
    if(length(celltype_tmp)>1){                          
      res_list_2 <- list()
      SNP_total_list <- list()
      for(c1 in 1:length(celltype_tmp)){
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

        ld_mat <- ieugwasr::ld_matrix(unique(harmonise_dat$SNP))
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
      
      ld_mat_SNP <- str_split_fixed(rownames(ld_mat),'_',2)[,1]
      harmonise_dat <- subset(harmonise_dat,SNP %in% ld_mat_SNP)
      harmonise_dat$celltype <- str_split_fixed(harmonise_dat$exposure,'_',2)[,2]
 
     
      beta_for_liml_PCA <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
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
      
      
      p1 <- nrow(bx)
      
      Phi1 <- ((rowSums(abs(bx))/sy)%*%t(rowSums(abs(bx))/sy))*ld
      r1 <- which(cumsum(prcomp(Phi1,scale=FALSE)$sdev^2/sum((prcomp(Phi1,scale=FALSE)$sdev^2)))>0.999)[1]
      lambda1 <- sqrt(p1)*prcomp(Phi1,scale=FALSE)$rotation[,1:r1]
      evec <- eigen((t(lambda1)%*%lambda1))$vectors
      eval <- eigen((t(lambda1)%*%lambda1))$values
      lambda1 <- lambda1%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
      dim(lambda1) <- c(p1,r1)
      
      
      BMA_dat_list <- list()
      for(mm1 in 1:length(celltype_tmp)){
        tmp_data_BMA1 <- pca.transform(beta = bx[,mm1],se = sx[,mm1],n = nx[1],ld = ld,by11 = beta_y,
                                      sy11 = sy,ny11 = ny,lambda11 = lambda1,r11=r1)
        BMA_dat_list[[mm1]] <- tmp_data_BMA1
      }
      BMA_by <- pca.transform(beta = beta_y,se = sy,n = ny,ld = ld,by11 = beta_y,
                              sy11 = sy,ny11 = ny,lambda11 = lambda1,r11=r1)
      exposure.id <- celltype_tmp
      pca.bx <- do.call(cbind,BMA_dat_list)
      colnames(pca.bx) <- exposure.id;
      pca.by <- BMA_by
      snps <- str_split_fixed(rownames(ld),'_',2)[,1]
      
      
      amd_nmr_input=new("mvMRInput", betaX = pca.bx, betaY = pca.by, snps=snps, exposure=exposure.id, outcome = "cad")
      BMA_output=summarymvMR_SSS(amd_nmr_input,kmin=1,kmax=10, prior_prob=0.1, max_iter=1000)
      
      BMA_pp <- function(kx){sort(BMA_output@pp,decreasing=TRUE)[1:length(celltype_tmp)][[kx]]} 
      BMA_pp <- sapply(1:length(celltype_tmp),BMA_pp)
      BMA_models <- function(kx1){exposure.id[as.numeric(strsplit(names(sort(BMA_output@pp,decreasing=TRUE)),split=",")[[kx1]])]} 
      sapply(1:length(celltype_tmp),BMA_models)
      round(BMA_pp,3)
      marg.inc <- cbind(exposure.id[order(BMA_output@pp_marginal, decreasing=TRUE)],round(BMA_output@pp_marginal[order(BMA_output@pp_marginal, decreasing=TRUE)],3))
      colnames(marg.inc) <- c("tissue","marg. inc. prob.")
    }
  }
}























