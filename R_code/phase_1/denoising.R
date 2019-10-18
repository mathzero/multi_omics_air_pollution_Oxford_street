rm(list=ls())
work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)


#############
library(lme4)
library(dplyr)
#############


#################################################################
#################################################################

###### mRNA #######

#load data set
mRNA <- readRDS("Project6/Data/phase_2/mRNA.RDS")

#Create empty denoised data set with colnames
mRNA.denoised <- mRNA
mRNA.denoised[,59:ncol(mRNA.denoised)] <- NA

#loop to denoise
for (i in 59:ncol(mRNA)){
  print(colnames(mRNA[i]))
  fml <- paste(colnames(mRNA[i]), "~ (1|array) + (1|batch) + Age + sex + BMI + Group")
  preds <- lmer(as.formula(fml), data = mRNA, REML = TRUE)
  mRNA.denoised[,i] <- resid(preds)
}

saveRDS(mRNA.denoised, "mRNA_denoised.rds")


#################################################################
#################################################################

###### miRNA #######

miRNA <- readRDS("Project6/Data/phase_2/miRNA.RDS")

miRNA.denoised <- miRNA
miRNA.denoised[,53:ncol(miRNA.denoised)] <- NA

#loop to denoise
for (i in 53:ncol(miRNA)){
  print(colnames(miRNA[i]))
  fml <- paste(colnames(miRNA[i]), "~ (1|Array) + (1|Batch) + Age + sex + BMI + Group")
  preds <- lmer(as.formula(fml), data = miRNA, REML = TRUE)
  miRNA.denoised[,i] <- resid(preds)
}

saveRDS(miRNA.denoised, "miRNA_denoised.rds")

#################################################################
#################################################################

###### metabolites #######

metabolites <- readRDS("Project6/Data/phase_2/metabolites.RDS")

metabolites.denoised <- metabolites
metabolites.denoised[,51:ncol(metabolites.denoised)] <- NA

#loop to denoise
for (i in 51:ncol(metabolites)){
  print(colnames(metabolites[i]))
  fml <- paste(colnames(metabolites[i]), "~ Age + sex + BMI + Group")
  preds <- lmer(as.formula(fml), data = metabolites, REML = TRUE)
  metabolites.denoised[,i] <- resid(preds)
}

saveRDS(metabolites.denoised, "metabolites_denoised.rds")

#################################################################
#################################################################

###### adductomics #######

adductomics <- readRDS("Project6/Data/phase_2/adductomics.RDS")

adductomics.denoised <- adductomics
adductomics.denoised[,61:ncol(adductomics.denoised)] <- NA

#loop to denoise
for (i in 61:ncol(adductomics)){
  print(colnames(adductomics[i]))
  fml <- paste(colnames(adductomics[i]), "~ Age + sex + BMI + Group")
  preds <- lmer(as.formula(fml), data = adductomics, REML = TRUE)
  adductomics.denoised[,i] <- resid(preds)
}

saveRDS(adductomics.denoised, "adductomics_denoised.rds")
