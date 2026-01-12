library(data.table)
library(vcfR)
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)




white_blood <- fread('UKB_dat1.csv')
identical(white_blood$eid,white_blood_2$eid)
total_data <- cbind(white_blood,white_blood_2[,-1])


diagnose_data <- fread('diagnose_data.csv')
attend_time <- read.csv('atend_time.csv')



diagnose_data_2 <- diagnose_data[match(rownames(ICD10_DM),diagnose_data$eid),]
identical(diagnose_data_2$eid,ICD10_DM$eid)
rownames(diagnose_data_2) <- diagnose_data_2$eid
identical(as.numeric(rownames(diagnose_data_2)),as.numeric(ICD10_DM$eid))
diagnose_data_2 <- diagnose_data_2[,-1]
diagnose_data_2 <- apply(diagnose_data_2,2,as.Date)
diagnose_data_2 <- data.frame(diagnose_data_2)

ethnic <- fread('ethnic.csv')
ethnic <- subset(ethnic,ethnic=='Euro')


total_data <- subset(total_data,eid %in% ethnic$eid)


load('DM_gp.Rdata')
gp_data <- data.frame(gp_data)
gp_data$attend_time <- attend_time$p53_i0[match(gp_data$eid,attend_time$eid)]
sum(is.na(gp_data$attend_time))
gp_data$time <- as.numeric(gp_data$time)
gp_data$check <- gp_data$attend_time-gp_data$time

con_dat <- gp_data[!(gp_data$eid %in% gp_data_com$eid),]




length(unique(gp_data_com$eid))
gp_data_com <- subset(gp_data_com,check<0)
total_data_test <- subset(total_data,eid %in% c(gp_data_com$eid,con_dat$eid))

tmp1 <- data.frame(table(gp_data$eid))
gp_data_uniqueID <- tmp1$Var1[tmp1$Freq<2]



total_data_2 <- total_data[total_data$eid %in% gp_data_com$eid,]
rownames(total_data_2) <- total_data_2$eid
name_id_1 <- c('p30000_i0',
               'p30120_i0','p30130_i0','p30140_i0',
               'p30150_i0','p30160_i0','p30180_i0','p30190_i0','p30200_i0',
               'p30210_i0','p30220_i0'
)

covar <- c('p30001_i0','p30002_i0','p30003_i0','p30004_i0',
           'p30121_i0','p30122_i0','p30123_i0','p30124_i0',
           'p30131_i0','p30132_i0','p30133_i0','p30134_i0',
           'p30141_i0','p30142_i0','p30143_i0','p30144_i0',
           'p30151_i0','p30152_i0','p30153_i0','p30154_i0',
           'p30161_i0','p30162_i0','p30163_i0','p30164_i0',
           'p30181_i0','p30182_i0','p30183_i0','p30184_i0',
           'p30191_i0','p30192_i0','p30193_i0','p30194_i0',
           'p30201_i0','p30202_i0','p30203_i0','p30204_i0',
           'p30211_i0','p30212_i0','p30213_i0','p30214_i0',
           'p30221_i0','p30222_i0','p30223_i0','p30224_i0'
)
total_data_counts <- total_data_2[,c(name_id_1)]
total_data_covar <- total_data_2[,c(covar)]
for(i in 1:11){
  total_data_covar[,4*(i-1)+1] <- as.numeric(total_data_covar[,4*(i-1)+1])
  total_data_covar[,4*(i-1)+2] <- as.numeric(as.Date(total_data_covar[,4*(i-1)+2]))
  total_data_covar[,4*(i-1)+3] <- as.factor(total_data_covar[,4*(i-1)+3])
  total_data_covar[,4*(i-1)+4] <- as.factor(total_data_covar[,4*(i-1)+4])
}




ageDF <- data.frame(fread('eid_BMI_Age_Sex.csv'))
ageDF_2 <- ageDF[match(rownames(total_data_covar),ageDF$eid),]

DF2 <- data.frame(ageDF_2[,c(2,3,4)],total_data_covar)
rownames(DF2) <- rownames(total_data_covar)
DF2 <- na.omit(DF2)
total_data_counts_2 <- total_data_counts[rownames(total_data_counts) %in% rownames(DF2),]
identical(rownames(total_data_counts_2),rownames(DF2))

counts_corr_list <- list()
num2 <- 1



for(j in 0:8){
  comp_ID <- j
  comp_Disease <- gp_data_com[gp_data_com$diseaseID %in% comp_ID,]
  data_DF1 <- rbind(con_dat,comp_Disease)
  data_DF1$diseasetype <- 0
  data_DF1$diseasetype[data_DF1$diseaseID %in% comp_ID] <- 1
  for(i in 1:11){
  DF3 <- data.frame(counts=total_data_counts_2[,i],DF2[,c(1,2,3,4*(i-1)+4,4*(i-1)+5,4*(i-1)+6,4*(i-1)+7)])
  DF4 <- DF3[match(data_DF1$eid,rownames(DF3)),]
  DF4 <- na.omit(DF4)
  data_DF2 <- data_DF1[match(rownames(DF4),data_DF1$eid),]
  data_DF3 <- data.frame(y=data_DF2$diseasetype,DF4)
  formula <- reformulate(termlabels = colnames(data_DF3)[-1], response = "y")
  model <- lm(formula,data=data_DF3)
  beta1<-summary(model)$coefficient[2,"Estimate"]
  se1<-summary(model)$coefficient[2,"Std. Error"]
  p1<-summary(model)$coefficient[2,"Pr(>|t|)"]
  lower1<-beta1-1.96*se1
  upper1<-beta1+1.96*se1
  OR1<-round(exp(beta1),3)
  ORlower1<-round(exp(lower1),3)
  ORupper1<-round(exp(upper1),3)
  counts_corr_list[[num2]] <- cbind(i,j,beta1,se1,p1,lower1,upper1,OR1,ORlower1,ORupper1)
  num2 <- num2+1
  }
}


corr_data <- do.call(rbind,counts_corr_list)
corr_data <- data.frame(corr_data)
corr_data$beta1[abs(corr_data$p1)>0.05] <- 0


corr_data2 <- corr_data[,c('i','j','beta1')]
corr_data3 <- dcast(corr_data2,formula = i~j)
corr_data3 <- corr_data3[,-1]
colnames(corr_data3) <- c('coma','ketoacidosis','renal complications',
                'ophthalmic complications','neurological complications',
                'peripheral circulatory complications','other specified complications',
                'multiple complications','unspecified complications')

library(uwot)
set.seed(1111)
tmp_data1 <- data.frame(umap(t(corr_data3),n_neighbors = 5))
rownames(tmp_data1) <- c('coma','ketoacidosis','renal complications',
                         'ophthalmic complications','neurological complications',
                         'peripheral circulatory complications','other specified complications',
                         'multiple complications','unspecified complications')
tmp_data2 <- tmp_data1[c(1,2,3,4,5,6),] 
tmp_data2$comp_name <- rownames(tmp_data2)
tmp_data2$name2 <- c('coma','Diabetic ketoacidosis','Diabetic nephropathy',
                     'Diabetic retinopathy',
                     'Diabetic neuropathy','Diabetic peripheral angiopathy')
UMAP_data <- tmp_data2

set.seed(1111)
library(cluster)
pam_result <- pam(tmp_data2[,1:2], k = 3)
pam_result$clustering









