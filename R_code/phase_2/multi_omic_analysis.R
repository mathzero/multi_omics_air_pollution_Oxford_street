rm(list=ls())
work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)

library(dplyr)

################################################
################################################

####  CREATING A COMBINED DATA SET OF ALL SIGNIFICANT OMICS ###

################################################
################################################

##### Importing lists of significant OMICS #####

significant_mRNA <- read.csv("Results/phase_1/results_mRNA/05_sensitivity/significant_mRNAs.csv")
significant_miRNA <- read.csv("Results/phase_1/results_miRNA/05_sensitivity/significant_miRNA_denoised.csv")
significant_metabolites <- read.csv("Results/phase_1/results_metabolomics/significant_metabolites.csv")
significant_adductomics <- read.csv("Results/phase_1/results_adductomics/significant_adductomics_denoised_cap05.csv")

##### Creating lists of significant OMICS from each data set

mRNA.sig <- as.character(significant_mRNA[,2])
miRNA.sig <- as.character(significant_miRNA[,2])
metab.sig <- as.character(significant_metabolites[,2])
add.sig <- as.character(significant_adductomics[,2])


################################################
################################################

#Read Covariates and exposures 
Full_covar_add <- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
                             sep = "\t", header = T)

### Impute means for missing ages ###

Full_covar_add$Age[is.na(Full_covar_add$Age)] <- mean(Full_covar_add$Age, na.rm = T)


for(i in 1:60){
  print(paste(colnames(Full_covar_add[i]),sum(is.na(Full_covar_add[,i]))))
}

######## Reading in data frame for metabolites #############


#Read imputed peaks
load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")
metabolites <- Full_peaks_imputed
colnames(metabolites)[2] <- "Group"
colnames(adductomics)[48] <- "id"

metabolites.sig.dat <- Full_peaks_imputed[,c("XP_ID", metab.sig)]




################################################

####  Reading in and transposing data frame for mRNA #######

OXF2_Tx_mRNA_procd_bkgrcorr_compld <- read.csv("Data/OXF2_Tx_mRNA_procd_bkgrcorr_compld.txt", sep="")

mRNA <- as.data.frame(t(OXF2_Tx_mRNA_procd_bkgrcorr_compld))
mRNA <- mRNA[2:nrow(mRNA),]
remove(OXF2_Tx_mRNA_procd_bkgrcorr_compld)


mRNA$idfull <- rownames(mRNA)
mRNA$id <- as.numeric(substring(mRNA$idfull,8,10))

mRNA$Time_Point <- as.factor(ifelse(grepl("_4", mRNA$idfull), "T2 (4PM)", 
                                    ifelse(grepl("_8", mRNA$idfull),"T1 (8AM)","T3 (24hrs)")))

numbers1 <- ifelse(grepl("Hy_8", mRNA$idfull), 11, 
                   ifelse(grepl("Hy_4", mRNA$idfull), 12, 
                          ifelse(grepl("Hy_24", mRNA$idfull), 13, 
                                 ifelse(grepl("Ox_8", mRNA$idfull), 21, 
                                        ifelse(grepl("Ox_4", mRNA$idfull), 22,23)))
                   )
)

mRNA$XP_ID <- as.numeric(paste(mRNA$id,numbers1, sep = ""))
mRNA <- mRNA[,c(30925:30928,1:30924)]



#### Cleaning mRNA tech confounders code #####

OXF2_Tx_mRNA_LabCode <- read.csv("Data/OXF2_Tx_mRNA_LabCode.csv")

mRNA.tech <- OXF2_Tx_mRNA_LabCode 

mRNA.tech$id <- as.numeric(substring(mRNA.tech$original.code,3,5))

numbers2 <- ifelse(grepl("1-8am", mRNA.tech$original.code), 11, 
                   ifelse(grepl("1-4pm", mRNA.tech$original.code), 12, 
                          ifelse(grepl("1-24h", mRNA.tech$original.code), 13, 
                                 ifelse(grepl("2-8am", mRNA.tech$original.code), 21, 
                                        ifelse(grepl("2-4pm", mRNA.tech$original.code), 22,23)))))

mRNA.tech$XP_ID <- as.numeric(paste(mRNA.tech$id,numbers2, sep = ""))

mRNA.tech <- mRNA.tech[,c(20:21,1:19)] #reordering variables
mRNA.tech$array <- as.integer(substring(mRNA.tech$Array.ID,1,5)) #extracting array

mRNA.tech <- mRNA.tech[,c(1:5,7:9,11:12,22)] #keeping only relevant variables
mRNA.tech <- mRNA.tech[mRNA.tech$Array.ID != "array not okay",] # removing dirty observations


##### binding mRNA and mRNA.tech ####

mRNA.full <-  inner_join(mRNA.tech, mRNA, by = "XP_ID")

#removing duplicate observation
mRNA <- mRNA.full[mRNA.full$idfull != "OXF2_RS026_Hy_24.1",]


####  indexing to select significant mRNAs #######

mRNA.sig.dat <- mRNA[,c("XP_ID", mRNA.sig)]


############################################################


######## Reading in data frame for miRNA #############

#Read miRNA
miRNA <- data.frame(t(read.table("Data/OXF2_Tx_miRNA_procd.txt")))
miRNA_target <- read.csv("Data/OXF2_Tx_miRNA_TargetFile.txt", sep="")

#changing "." in rownames to "-"
rownames <- rownames(miRNA)
rownames <- gsub("[.]", "-", rownames)

rownames(miRNA) <- rownames

miRNA$SampleCode <- rownames(miRNA) 

miRNA <- inner_join(miRNA_target, miRNA, by = "SampleCode")

miRNA$ID <- as.numeric(substring(miRNA$SampleCode,3,5))
miRNA$Time_Point <- as.factor(ifelse(grepl("4pm", miRNA$SampleCode), "T2 (4PM)", 
                                     ifelse(grepl("8am", miRNA$SampleCode),"T1 (8AM)","T3 (24hrs)")))
#hyde park = 1, Oxford =2
# treatment 1-3 = HP, treatment 4-6 = Oxford

numbers <- ifelse(miRNA$Treatment == 1, 11, 
                  ifelse(miRNA$Treatment == 2, 12, 
                         ifelse(miRNA$Treatment == 3, 13, 
                                ifelse(miRNA$Treatment == 4, 21, 
                                       ifelse(miRNA$Treatment == 5, 22,23)))
                  )
)

miRNA$XP_ID <- as.numeric(paste(miRNA$ID,numbers, sep = ""))

miRNA <- miRNA[,c(373:375,1:372)] 
####  indexing to select significant miRNAs #######

miRNA.sig.dat <- miRNA[,c("XP_ID", miRNA.sig)]

############################################################
######## Reading in data frame for adductomics #############

#Read adductomics
adductomics <- read.csv("Data/Adductomics_Oxford_Street_all_imputed.csv")
# removing batch 56 from adductomics as it apears to be duplication
adductomics <- adductomics[adductomics$batch != 56,]
rownames(adductomics) <- adductomics$XP_ID

#log transforming adductomics
adductomics[,18:ncol(adductomics)] <- log(adductomics[,18:ncol(adductomics)])


adductomics.sig.dat <- adductomics[,c("XP_ID",add.sig)]






########################################################################################################################
########################################################################################################################
########################################################################################################################



#### Joining all data sets #####
library(dplyr)

data.sig <- inner_join(miRNA.sig.dat,mRNA.sig.dat, by = "XP_ID")
data.sig <- inner_join(data.sig,metabolites.sig.dat, by = "XP_ID")
data.sig <- inner_join(data.sig,adductomics.sig.dat, by = "XP_ID")
data.sig <- inner_join(data.sig,Full_covar_add, by = "XP_ID")

sig.ids <- data.sig$XP_ID


# Data notes
 
# miRNA: Cols 2:7
# mRNA: Cols 8:283
# metabs: Cols 284:431
# adducts: Cols 432:434

write.csv(data.sig, "multiomic_significant_vars.csv")



########################################################################################################################
########################################################################################################################
########################################################################################################################

### Creating data sets for denoising ###

#adductomics
adductomics <- adductomics[match(sig.ids,adductomics$XP_ID),]
adductomics <- inner_join(Full_covar_add, adductomics, by = "XP_ID")
saveRDS(adductomics, "adductomics.RDS")

#miRNA
miRNA <- miRNA[match(sig.ids,miRNA$XP_ID),]
miRNA <- inner_join(Full_covar_add, miRNA, by = "XP_ID")
saveRDS(miRNA, "miRNA.RDS")

#metabolites
metabolites <- metabolites[match(sig.ids,metabolites$XP_ID),]
metabolites <- inner_join(Full_covar_add,metabolites, by = "XP_ID")
saveRDS(metabolites, "metabolites.RDS")

#mRNA
mRNA <- mRNA[match(sig.ids,mRNA$XP_ID),]
mRNA <- inner_join(Full_covar_add,mRNA, by = "XP_ID")
saveRDS(mRNA, "mRNA.RDS")


metabolites <- metabolites[,1:(ncol(metabolites)-2)]

metabolites[1:5,(ncol(metabolites)-2):ncol(metabolites)]










########################################################################################################################
########################################################################################################################
########################################################################################################################

#Experimenting with heatmaps

data.sig[,2:(160)] <- as.data.frame(scale(data.sig[,2:(160)]))

mat <- as.matrix(data.sig[data.sig$Per == "after_1",2:(160)])
cormat <- cor(mat)
heatmap(mat, Rowv=FALSE)



data.sig[order(Loc),]

library(gplots)
library(RColorBrewer)

hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)

#Visualising expression levels




png(file = "expression_levels.png", bg = "transparent", width = 2000, height = 2000, res =300)

heatmap.2(as.matrix(mat),
          density.info="none", 
          trace="none", 
          Colv=FALSE, 
          Rowv=FALSE, 
          dendrogram="none", 
          margin=c(10,10),
          keysize = 1,
          cexRow=0.3,
          cexCol=0.3,
          col = rev(hmcol),
          # now the separations:
          colsep=c(6,8,157))
dev.off()  # 318 bytes file in current directory

#heatmap
png(file = "omics_heatmap.png", bg = "transparent", width = 2500, height = 2500, res =300)
heatmap.2(as.matrix(cormat),
          density.info="none", 
          trace="none", 
          Colv=FALSE, 
          Rowv=FALSE, 
          dendrogram="none", 
          margin=c(10,10),
          keysize = 1,
          cexRow=0.3,
          cexCol=0.3,
          col = rev(hmcol),
          # now the separations:
          colsep=c(6,8,157),
          rowsep=c(6,8,157))
dev.off()  # 318 bytes file in current directory


dat.after1 <- data.sig[data.sig$Per == "after_1",]
dat.after1[order(dat.after1$Loc),"Loc"]

### Expression map separated by location
png(file = "omics_expression_by_location.png", bg = "transparent", width = 2500, height = 2500, res =300)
heatmap.2(as.matrix(dat.after1[order(dat.after1$Loc),2:160]),
          density.info="none", 
          trace="none", 
          Colv=FALSE, 
          Rowv=FALSE, 
          dendrogram="none", 
          margin=c(10,10),
          keysize = 1,
          cexRow=0.3,
          cexCol=0.3,
          col = rev(hmcol),
          xlab  = "Biomarker",
          ylab = "Observation, ordered by location: OXFORD STREET / HYDE PARK",
          # now the separations:
          colsep=c(6,8,157),
          rowsep=39)

dev.off() 

table(dat.after1$Group)

### Expression map separated by health status
png(file = "omics_expression_by_health_status.png", bg = "transparent", width = 2500, height = 2500, res =300)
heatmap.2(as.matrix(dat.after1[order(dat.after1$Group),2:160]),
          density.info="none", 
          trace="none", 
          Colv=FALSE, 
          Rowv=FALSE, 
          dendrogram="none", 
          margin=c(10,10),
          keysize = 1,
          cexRow=0.3,
          cexCol=0.3,
          col = rev(hmcol),
          xlab  = "Biomarker",
          ylab = "Participant, ordered by health status: IHD / HEALTHY / COPD",
          # now the separations:
          colsep=c(6,8,157),
          rowsep=c(32,58))
dev.off() 

?heatmap.2

########################################################################################################################
########################################################################################################################
########################################################################################################################



























# Univariate analysis of significant miRNAs against all metabolites


#Load miRNA
miRNA.1 = scale(data.sig[,unlist(miRNA.sig)])

load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")

#Load metabolites
metabolites <- as.data.frame(Full_peaks_imputed[,c("XP_ID", metab.sig)])
metabolites[,2:ncol(metabolites)] <- scale(metabolites[,2:ncol(metabolites)])
rownames(metabolites) <- metabolites$XP_ID
metabolites <- as.data.frame(metabolites[sig.ids,2:ncol(metabolites)])
Beta = Pvalues = matrix(NA, nrow = ncol(miRNA.1),ncol = ncol(metabolites))
                        

ncol(miRNA.1) * ncol(metabolites) # number of tests
t0 = Sys.time()
for (i in 1:ncol(miRNA.1)) {
  print(i)
  for (j in 1:ncol(metabolites)) {
    model1 = lm(miRNA.1[, i] ~ metabolites[, j])
    Beta[i, j] = coefficients(model1)["metabolites[, j]"]
    Pvalues[i, j] = summary(model1)$coefficients["metabolites[, j]",
                                                 "Pr(>|t|)"]
  }
}
t1 = Sys.time()
print(t1 - t0)
rownames(Pvalues) = rownames(Beta) = colnames(miRNA.1)
colnames(Pvalues) = colnames(Beta) = colnames(metabolites)


#Save results
dir.create("Results_univ", showWarnings = FALSE)
saveRDS(Pvalues, "Results_univ/pvalues_miRNA-metab.rds")
saveRDS(Beta, "Results_univ/betas_miRNA-metab.rds")

plot(Pvalues)
sum(p.adjust(as.vector(Pvalues), method = "BH") <
      0.05)
