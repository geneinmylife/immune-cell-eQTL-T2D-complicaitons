#steiger

setwd("/Users/dingyilan/Documents/dm")
dat <- read.csv("dat_dm.csv")
dat$samplesize.exposure <- 91
dat$samplesize.outcome <-dat$ncase+dat$ncontrol
dat$rsq.exposure <- (get_r_from_pn(p=dat$pval.exposure,n=dat$samplesize.exposure))^2
dat$rsq.outcome <- (get_r_from_lor(
  lor=dat$beta.outcome,
  af=dat$eaf.outcome,
  ncase=dat$ncase,
  ncontrol=dat$ncontrol,
  prevalence=dat$ncase/dat$samplesize.outcome))^2
st <- psych::r.test( 
  n = dat$samplesize.exposure, 
  n2 = dat$samplesize.outcome, 
  r12 = sqrt(dat$rsq.exposure), 
  r34 = sqrt(dat$rsq.outcome))
dat$steiger_dir <- dat$rsq.exposure >dat$rsq.outcome
dat$steiger_pval <- pnorm(-abs(st$z)) * 2 

write.csv(dat,'dm_steiger.csv')
