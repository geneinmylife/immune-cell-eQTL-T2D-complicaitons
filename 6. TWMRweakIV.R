library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)



#TWMR


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



TWMR_res_list <- list()
num1 <- 1
file_list <- list()
for(d1 in 1:length(disease_info$name)){
  tmp_res_d1 <- subset(res_data_total,outcome==disease_info$name[d1])
  tmp_gene_list_1 <- unique(tmp_res_d1$id)
  R9_data <- R9_data_list_total[[d1]]
  for(g1 in 1:length(tmp_gene_list_1)){
    tmp_res_d2 <- subset(tmp_res_d1,id==tmp_gene_list_1[g1])
    celltype_tmp <- unique(tmp_res_d2$celltype)
    if(length(celltype_tmp)>1){                           
      SNP_tmp_x1 <- read.table(paste0('SNP_pruned/',disease_info$name[d1],'_',tmp_gene_list_1[g1],'.prune.in'),
                               header = F)
      SNP_ID_tmp1 <- SNP_tmp_x1$V1
      
      
      #match Beta SE
      SNP_total_list_2 <- list()
      for(c1 in 1:length(celltype_tmp)){
        tmp_gene_data <- fread(paste0('MV_IVW_PCA/SNP_celltype/',celltype_tmp[c1],'_data_for_MV_IVW_PCA.csv'))
        df_separated_x1 <- tmp_gene_data
        SNP_total_list_2[[c1]] <- df_separated_x1
      }
      SNP_total_dat_2 <- do.call(rbind,SNP_total_list_2)
      df_separated_tmp_2 <- subset(SNP_total_dat_2,Gene==tmp_gene_list_1[g1])
      df_separated_tmp_2 <- subset(df_separated_tmp_2,ID %in% SNP_ID_tmp1)
      
      df_separated_tmp_2$unique_SNPID <- paste(df_separated_tmp_2$ID,df_separated_tmp_2$REF,df_separated_tmp_2$ALT,sep = '_')
      df_separated_tmp_2$unique_SNPID2 <- paste(df_separated_tmp_2$ID,df_separated_tmp_2$ALT,df_separated_tmp_2$REF,sep = '_')
      df_separated_tmp_2 <- subset(df_separated_tmp_2,unique_SNPID %in% R9_data$unique_SNPID | unique_SNPID %in% R9_data$unique_SNPID2)
      
      
      
      SNP_ID_use <- unique(df_separated_tmp_2$ID)
      
      #calculate ld matrix
      ld_mat <- ieugwasr::ld_matrix(SNP_ID_use,plink_bin = 'software/plink/plink.exe',
                                    bfile = 'data/g1000_eur/g1000_eur')
      
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
      
      beta_for_TWMR <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      se_for_TWMR <- data.frame(matrix(0,nrow = nrow(ld_mat),ncol = length(celltype_tmp)))
      for(j in 1:length(celltype_tmp)){
        df_separated_tmp_3 <- data.frame(subset(df_separated_tmp_2,celltype==celltype_tmp[j]))
        rownames(df_separated_tmp_3) <- df_separated_tmp_3$SNPID_same_with_ld_mat
        df_separated_tmp_3 <- df_separated_tmp_3[rownames(ld_mat),]
        beta_for_TWMR[,j] <- df_separated_tmp_3$Beta2
        se_for_TWMR[,j] <- df_separated_tmp_3$SE
      }
      
      se_y <- R9_dat_tmp_1$sebeta
      beta_y <- R9_dat_tmp_1$Beta2
      
      
      #TWMR
      if(ncol(beta_for_TWMR)<=nrow(beta_for_TWMR)){
        Ngwas<-disease_info$ncase[d1]+disease_info$ncontrol[d1]
        N_eQTLs<-91
        out<-c("gene","alpha","SE","P","Nsnps","Ncelltype","unique")
        
        beta<-beta_for_TWMR
        
        beta <- as.matrix(beta)
        gamma <- as.matrix(beta_y)
        
        
        C<-ld_mat
        C<-as.matrix(C)
        
        S<-t(beta)%*%solve(C)%*%beta
        H<-(1-1/sqrt(3781))*S+(1/sqrt(3781))*diag(length(S[,1]))
        alpha<-solve(H)%*%(t(beta)%*%solve(C)%*%gamma)
        
        alpha<-as.vector(alpha)
        
        C_inv <- solve(C)
        GCG_inv <- t(beta) %*% solve(C) %*% beta
        GCG_inv<-(1-1/sqrt(3781))*GCG_inv+(1/sqrt(3781))*diag(length(GCG_inv[,1]))
        GCG_inv<-solve(GCG_inv)
        
        
        df_dg <- GCG_inv %*% t(beta) %*% C_inv
        df_dG <- (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))))) + ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
        J <- cbind(df_dG, df_dg)
        
        SEs<-c(rep(1/sqrt(N_eQTLs),length(beta[1,])*length(beta[,1])),rep(1/sqrt(Ngwas),length(gamma[,1])))
        R<-diag(length(beta[1,])+1)
        Sigma <- (SEs %*% t(SEs)) * (C %x% R)   
        V <- J %*% Sigma %*% t(J)
        
        
        out_list <- list()
        for(k1 in 1:ncol(beta)){
          se<- sqrt(V[k1,k1])
          N=nrow(beta)
          Ngene=ncol(beta)
          Z<-alpha[k1]/se
          pval<-2*pnorm(abs(Z),lower.tail=FALSE)
          line<-c(celltype_tmp[k1],alpha[k1],se,pval,N,Ngene,paste(disease_info$name[d1],celltype_tmp[k1],tmp_gene_list_1[g1],sep = '_'))
          out_list[[k1]] <- line
        }
        out_total <- do.call(rbind,out_list)

        TWMR_res_list[[num1]] <- out_total
        num1 <- num1+1
      }
    }
  }
}






