library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(TwoSampleMR)
library(ieugwasr)
library(MVMR)
library(mvmr.weakiv)

#Weak instrument robust two-sample MVMR


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

MVMR_res_list <- list()

num1 <- 1
num2 <- 1
for(d1 in 1:length(disease_info$name)){#
  tmp_res_d1 <- subset(res_data_total,outcome==disease_info$name[d1])
  tmp_gene_list_1 <- unique(tmp_res_d1$id)
  R9_data <- R9_data_list_total[[d1]]
  for(g1 in 1:length(tmp_gene_list_1)){#
    tmp_res_d2 <- subset(tmp_res_d1,id==tmp_gene_list_1[g1])
    celltype_tmp <- unique(tmp_res_d2$celltype)
    if(length(celltype_tmp)>1){                           #only pleiotropy
      res_list_2 <- list()
      SNP_total_list <- list()
      rho_tmp <- rho[celltype_tmp,celltype_tmp]
      
      
      #Extract all independent IVs after clumping
      for(c1 in 1:length(celltype_tmp)){
        tmp_gene_data <- fread(paste0(celltype_tmp[c1],'_clumped_steiger_data.csv'))
        df_separated <- tmp_gene_data
        SNP_total_list[[c1]] <- df_separated
      }
      SNP_total_dat <- do.call(rbind,SNP_total_list)
      SNP_total_dat$disease <- str_split_fixed(SNP_total_dat$unique,'_ENSG',2)[,1]
      SNP_total_dat$disease <- gsub('.*_','',SNP_total_dat$disease)
      SNP_total_dat <- subset(SNP_total_dat,disease==disease_info$name[d1])
      SNP_total_dat$Gene <- str_split_fixed(SNP_total_dat$exposure,'_',2)[,1]
      df_separated_tmp <- subset(SNP_total_dat,Gene==tmp_gene_list_1[g1])
      df_separated_tmp$unique_SNPID <- paste(df_separated_tmp$rsid,df_separated_tmp$other_allele.exposure,df_separated_tmp$effect_allele.exposure,sep = '_')
      df_separated_tmp$unique_SNPID2 <- paste(df_separated_tmp$rsid,df_separated_tmp$effect_allele.exposure,df_separated_tmp$other_allele.exposure,sep = '_')
      #Keep SNPs present in both FinnGen and DICE
      df_separated_tmp <- subset(df_separated_tmp,unique_SNPID %in% R9_data$unique_SNPID | unique_SNPID %in% R9_data$unique_SNPID2)
      SNP_ID_tmp1 <- unique(df_separated_tmp$rsid)
      
      
      
      #match Beta SE
      SNP_total_list_2 <- list()
      for(c1 in 1:length(celltype_tmp)){
        tmp_gene_data <- fread(paste0(celltype_tmp[c1],'_data_for_MV_IVW_PCA.csv'))
        df_separated_x1 <- tmp_gene_data
        SNP_total_list_2[[c1]] <- df_separated_x1
      }
      SNP_total_dat_2 <- do.call(rbind,SNP_total_list_2)
      df_separated_tmp_2 <- subset(SNP_total_dat_2,Gene==tmp_gene_list_1[g1])
      df_separated_tmp_2 <- subset(df_separated_tmp_2,ID %in% SNP_ID_tmp1)
      tmp1 <- data.frame(table(df_separated_tmp_2$ID))
      tmp2 <- subset(tmp1,Freq>=length(celltype_tmp))
      df_separated_tmp_2 <- subset(df_separated_tmp_2,ID %in% tmp2$Var1)
      
      df_separated_tmp_2$unique_SNPID <- paste(df_separated_tmp_2$ID,df_separated_tmp_2$REF,df_separated_tmp_2$ALT,sep = '_')
      df_separated_tmp_2$unique_SNPID2 <- paste(df_separated_tmp_2$ID,df_separated_tmp_2$ALT,df_separated_tmp_2$REF,sep = '_')
      
      df_separated_tmp_2 <- subset(df_separated_tmp_2,unique_SNPID %in% R9_data$unique_SNPID | unique_SNPID %in% R9_data$unique_SNPID2)
      
      
      
      SNP_ID_use <- unique(df_separated_tmp_2$ID)
      
      #calculate ld matrix
      ld_mat <- ieugwasr::ld_matrix(SNP_ID_use,plink_bin = 'software/plink/plink.exe',
                                    bfile = 'data/g1000_eur/g1000_eur')
      remove1 <- c()
      for(ld1 in 1:nrow(ld_mat)){
        tmp_ld1 <- ld_mat[,ld1]
        for(ld2 in ld1:length(tmp_ld1)){
          if(tmp_ld1[ld2]==1 & ld1 != ld2){
            remove1 <- c(remove1,ld2)
          }
        }
      }
      if(length(remove1)>0){
        ld_mat <- ld_mat[-remove1,-remove1]
      }
      ld_mat_SNP <- str_split_fixed(rownames(ld_mat),'_',2)[,1]
      
      df_separated_tmp_2 <- subset(df_separated_tmp_2,ID %in% ld_mat_SNP)
      
      df_separated_tmp_2$SNPID_same_with_ld_mat <- df_separated_tmp_2$unique_SNPID
      df_separated_tmp_2$Beta2 <- df_separated_tmp_2$Beta
      df_separated_tmp_2$Beta2[df_separated_tmp_2$unique_SNPID2 %in% rownames(ld_mat)] <- -df_separated_tmp_2$Beta2[df_separated_tmp_2$unique_SNPID2 %in% rownames(ld_mat)] 
      df_separated_tmp_2$SNPID_same_with_ld_mat[df_separated_tmp_2$unique_SNPID2 %in% rownames(ld_mat)] <- df_separated_tmp_2$unique_SNPID2[df_separated_tmp_2$unique_SNPID2 %in% rownames(ld_mat)]
      
      R9_dat_tmp_1 <- subset(R9_data,unique_SNPID %in% df_separated_tmp_2$SNPID_same_with_ld_mat | 
                               unique_SNPID2 %in% df_separated_tmp_2$SNPID_same_with_ld_mat)
      R9_dat_tmp_1$SNPID_same_with_ld_mat <- R9_dat_tmp_1$unique_SNPID
      R9_dat_tmp_1$Beta2 <- R9_dat_tmp_1$beta
      R9_dat_tmp_1$Beta2[R9_dat_tmp_1$unique_SNPID2 %in% rownames(ld_mat)] <- -R9_dat_tmp_1$Beta2[R9_dat_tmp_1$unique_SNPID2 %in% rownames(ld_mat)]
      R9_dat_tmp_1$SNPID_same_with_ld_mat[R9_dat_tmp_1$unique_SNPID2 %in% rownames(ld_mat)] <- R9_dat_tmp_1$unique_SNPID2[R9_dat_tmp_1$unique_SNPID2 %in% rownames(ld_mat)]
      rownames(R9_dat_tmp_1) <- R9_dat_tmp_1$SNPID_same_with_ld_mat
      R9_dat_tmp_1 <- R9_dat_tmp_1[rownames(ld_mat),]
      identical(rownames(R9_dat_tmp_1),rownames(ld_mat))
      
      beta_for_MVMR <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      se_for_MVMR <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      for(j in 1:length(celltype_tmp)){
        df_separated_tmp_3 <- data.frame(subset(df_separated_tmp_2,celltype==celltype_tmp[j]))
        rownames(df_separated_tmp_3) <- df_separated_tmp_3$SNPID_same_with_ld_mat
        df_separated_tmp_3 <- df_separated_tmp_3[rownames(ld_mat),]
        beta_for_MVMR[,j] <- df_separated_tmp_3$Beta2
        se_for_MVMR[,j] <- df_separated_tmp_3$SE
      }
      
      se_y <- R9_dat_tmp_1$sebeta
      beta_y <- R9_dat_tmp_1$Beta2
      
      #mvmr.weakiv
      if(ncol(beta_for_MVMR)<=nrow(beta_for_MVMR)){
        f <- try(res <- mvmr.weakiv(bx=beta_for_MVMR,
                                    by=beta_y,
                                    sx=se_for_MVMR,
                                    sy=se_y,
                                    ld=ld_mat,
                                    nx=91,
                                    ny=disease_info$ncase[d1]+disease_info$ncontrol[d1],
                                    cor.x=rho_tmp,
                                    max.search=0.3,len.search=30),silent = TRUE)
        
        if(class(f)=='try-error'){
          next
        } else{
          MVMR_res_list[[num1]] <- res
          num1 <- num1+1
        }
      }
    }
  }
}

