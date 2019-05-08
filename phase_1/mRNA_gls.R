rm(list=ls())
work_dir <- "/rds/general/user/mw418/home/TDS_Project"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"


setwd(work_dir)

#load libraries
library(mvtnorm)
library(nlme)
library(doBy)
library(dplyr)
library(lme4)
library(parallel)


####  Reading in and transposing data frame #######

OXF2_Tx_mRNA_procd_bkgrcorr_compld <- read.csv("Data/OXF2_Tx_mRNA_procd_bkgrcorr_compld.txt", sep="")

mRNA <- as.data.frame(t(OXF2_Tx_mRNA_procd_bkgrcorr_compld))
mRNA <- mRNA[2:nrow(mRNA),]

##### Manipulating data to extract XP_ID ######

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

###############################################################
###############################################################


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
mRNA.full <- mRNA.full[mRNA.full$idfull != "OXF2_RS026_Hy_24.1",]

###### ###### ###### ######




###############################################################
###############################################################



#Read Covariates and exposures 
Full_covar_add <- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
                             sep = "\t", header = T)

#Read imputed peaks
load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")

# Merge information from the two datasets
# Make XP_ID the rownames
rownames(Full_covar_add) <- Full_covar_add$XP_ID
rownames(mRNA.full) <- mRNA.full$XP_ID
rownames(Full_peaks_imputed) <- Full_peaks_imputed$XP_ID


# select people who have all information (50)
select.idx <- intersect(rownames(Full_covar_add),rownames(mRNA.full))

# just keep these individuals (50 individuals, 3 measurements each)
Full_covar_add <- Full_covar_add[select.idx,]
mRNA.full <- mRNA.full[select.idx,]
Full_peaks_imputed <- Full_peaks_imputed[select.idx,]


# store variables of interest
data <- NULL
data <- Full_covar_add[,c("ID","XP_ID","Group","Loc","Per","Age", 
                          "sex", "BMI", "Group", "Caffeine")]
data <- cbind(data,Full_peaks_imputed[,c("Site","Time_Point")])
data <- left_join(data,mRNA.tech[,c("XP_ID","array","batch")], by = "XP_ID")
data$wave <- factor(as.numeric(factor(10*as.numeric(factor(data$Site))+
                                        as.numeric(factor(data$Time_Point)))))




# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")

# Omics 
ListOmics4 <- colnames(mRNA.full[16:ncol(mRNA.full)])

# limits (repeat with FALSE, 0.05, 0.1, 0.2)
capOMIC <- TRUE
capValOMIC <- 0.05


data$Age <- as.numeric(data$Age)
data$Caffeine <- as.factor(data$Caffeine)
data$batch <- as.integer(data$batch)


#Recode Group variable so that healthy becomes the reference
data$group2[data$Group == "HEALTHY"] <- "1HEALTHY"
data$group2[data$Group == "COPD"] <- "2COPD"
data$group2[data$Group == "IHD"] <- "3IHD"




###############################################################
###############################################################

######## FORMULA #########



NOmics <- length(ListOmics4)

fun.loop <- function(ListExpo2){
  for(MyExp in ListExpo2[1:5]){
    FileName <- paste("Results_",MyExp,sep="")
    t00 <- Sys.time()
    for(CurrOmic in 1:NOmics){
      
      # create exposure variables in the dataset
      # two variables (1 annual)
      eval(parse(text=paste("data$Expo <- Full_covar_add$",MyExp, "_o",sep=''))) 
      eval(parse(text=paste("data$Expo_annual <- Full_covar_add$",
                            MyExp,"_annual",sep='')))
      
      # set time 
      t0 <- Sys.time()
      # create omic variable in the dataset
      MyOmic <- ListOmics4[CurrOmic]
      print(paste("OMICS signal ", CurrOmic,"/",length(ListOmics4)," - ", MyOmic,sep=''))
      eval(parse(text=paste("data$Omic <- mRNA.full$",MyOmic,sep='')))
      
      #standardise OMIC measurement
      data$Omic <- data$Omic/sd(data$Omic)
      
      # correct extreme values
      # select values outside quantile region and set them to the extreme 
      # repeat for upper and lower quantile
      if(capOMIC){
        data$Omic[data$Omic>quantile(data$Omic,(1-capValOMIC))]=quantile(data$Omic,(1-capValOMIC)) # If omic value is in the >80th percentile, reset the value to the value of the 80th percentile
        data$Omic[data$Omic<quantile(data$Omic,capValOMIC)]=quantile(data$Omic,capValOMIC) # If omic value is in the <82th percentile, reset the value to the value of the 20th percentile
      }
      
      # create variable for deviation 
      # Use delta=0 if T1 (time before)
      # Use delta= backpack - annual otherwise (lagging effect)
      
      data$Expo.delta<-ifelse(data$Per=="before", 0, data$Expo-data$Expo_annual)
      
      # gls model for omic measurement
      # adjusted for age, sex, bmi, caffeine and group
      # multiple measurements per individual (1|ID)
      # variance covariance across wave (1|wave)
      # include just data points with no missing information
      
      preds <- lmer(Omic ~(1|array) + (1|batch), data = data, REML = TRUE)
      data$Omic <- resid(preds)
      
      M6 <- gls(I(Omic) ~ Expo_annual + Expo.delta + Age + sex + BMI + 
                  group2 + Caffeine,
                data=data[complete.cases(data),],
                correlation=corSymm(form=~1|ID),
                weights=varIdent(form=~1|wave),
                na.action=na.pass,
                control = list(singular.ok = TRUE))
      
      
      
      
      #Storing Results: coeffs and p-values separately
      if(MyOmic==ListOmics4[1]){
        
        resM6Pval<- summary(M6)$tTable[,"p-value"]
        resM6Coeff<- summary(M6)$tTable[,"Value"] 
        
        
      }else{
        
        
        resM6Pval<- rbind(resM6Pval,summary(M6)$tTable[,"p-value"])
        resM6Coeff<- rbind(resM6Coeff,summary(M6)$tTable[,"Value"])
        
      }
      
      
      # elapsed time
      t1 <- Sys.time()
      # user message
      print(paste("Time elapsed ",round(as.numeric(t1-t0,units="secs"),digits=2),
                  " secs", " -- Total : ", round(as.numeric(t1-t00,units="secs"),
                                                 digits=2),
                  ' secs - ',round(as.numeric(t1-t00,units="secs")/60,digits=2), 
                  ' mins', sep=''))
      #Reset the OMIC variable at the end of the loop
      data <- data[,!(names(data) %in% "Omic")]
      #Reset the exposure variable at the end of the loop
      data <- data[,!(names(data) %in% c("Expo_annual","Expo","Expo.delta"))]
    }
    
    # total time
    t2 <- Sys.time()
    # user message
    print(paste("Total time: ", round(as.numeric(t2-t00,units="secs"),digits=2),
                " secs - ", round(as.numeric(t2-t00,units="secs")/60,digits=2),
                " minutes - ", round(as.numeric(t2-t00,units="secs")/3600,digits=2),
                " hours",sep=''))
    
    # store results together
    mRNA_resM6Pval <- cbind(ListOmics4[1:NOmics],resM6Pval)
    colnames(mRNA_resM6Pval)[1] <- "OMIC_ID"
    mRNA_resM6Coeff <- cbind(ListOmics4[1:NOmics],resM6Coeff)
    colnames(mRNA_resM6Coeff)[1] <- "OMIC_ID"
    # save results 
    write.table(mRNA_resM6Pval, paste(FileName,"_mRNA_M6_Pvals.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    write.table(mRNA_resM6Coeff, paste(FileName, "_mRNA_M6_Coeff.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    
    
    
  }
}


##### Executing loop ####


# fun.loop(ListExpo2)




######  Parrallelised (uncomment to execute)  ######


# t0=Sys.time()
# no_cores=detectCores()
# cl <- makeCluster(no_cores, type="FORK")
# 
# 
# system.time(parLapply(cl, ListExpo2, fun.loop))
# 
# stopCluster(cl)
# t1=Sys.time()
# print(t1-t0)



#######################

### Multi-core parallelised

##### Executing loop ####


ids=as.character(cut(1:length(ListOMICs2), breaks = nchunks, labels = 1:nchunks))

t0=Sys.time()
no_cores=detectCores()
cl <- makeCluster(no_cores) 
clusterExport(cl, c("work_dir", "data", "Full_covar_add", "ListExpo2", "ListOMICs2", "mRNA"))
clusterEvalQ(cl, library(mvtnorm),library(nlme),library(doBy))

pvalues=parSapply(cl=cl, 1:nchunks, FUN=function(k){
  X_chunk=ListOMICs2[,ids==k]
  return(apply(X_chunk, 2, FUN = fun.loop))
})

stopCluster(cl)
t1=Sys.time()
print(t1-t0)






####Visualisations ####


work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/Results/phase_1/results_mRNA/05_sensitivity"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"


setwd(work_dir)




library(ggplot2)
library(ggrepel)





mRNA_NO2_imp <-cbind(read.table("Results_NO2_mRNA_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_NO2_mRNA_M6_Coeff.txt", header = TRUE)[4])
mRNA_PM10_imp <-cbind(read.table("Results_PM10_mRNA_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM10_mRNA_M6_Coeff.txt", header = TRUE)[4])
mRNA_PM25_imp <-cbind(read.table("Results_PM25_mRNA_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM25_mRNA_M6_Coeff.txt", header = TRUE)[4])
mRNA_PCNT_imp <-cbind(read.table("Results_PCNT_mRNA_M6_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_PCNT_mRNA_M6_Coeff.txt", header = TRUE)[3])
mRNA_CBLK_imp <-cbind(read.table("Results_CBLK_mRNA_M6_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_CBLK_mRNA_M6_Coeff.txt", header = TRUE)[3])


colnames(mRNA_NO2_imp) <- c("OMIC.ID","pval","beta")
colnames(mRNA_PM10_imp) <- c("OMIC.ID","pval","beta")
colnames(mRNA_PM25_imp) <- c("OMIC.ID","pval","beta")
colnames(mRNA_PCNT_imp) <- c("OMIC.ID","pval","beta")
colnames(mRNA_CBLK_imp) <- c("OMIC.ID","pval","beta")

# read results for each mRNAcription exposure# 
exposures.mRNA <- c("mRNA_NO2_imp","mRNA_PM10_imp","mRNA_PM25_imp","mRNA_PCNT_imp","mRNA_CBLK_imp")

# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")

mRNA.list <- list(mRNA_NO2_imp,mRNA_PM10_imp,mRNA_PM25_imp,mRNA_PCNT_imp,mRNA_CBLK_imp)
names(mRNA.list) <- exposures.mRNA

##############################################################################################
##############################################################################################

### Creating full data frame pf pvals for manhattan plots ####

mRNA.full.sig <- cbind(mRNA.list[[1]][1:2], mRNA.list[[2]][2],mRNA.list[[3]][2], mRNA.list[[4]][2], mRNA.list[[5]][2])

colnames(mRNA.full.sig) <- c("OMICid","NO2_pval","PM10_pval","PM25_pval","PCNT_pval","CBLK_pval")
write.csv(mRNA.full.sig, "mRNA_all_pvals.csv")

##############################################################################################
##############################################################################################



for (i in 1:length(exposures.mRNA)){
  mRNA.list[[i]]$pvalBH <- round(p.adjust(mRNA.list[[i]]$pval, method = "BH"),3)
  mRNA.list[[i]]$pvalFDR <- round(p.adjust(mRNA.list[[i]]$pval, method = "fdr"),3)
  print(exposures.mRNA[i])
  print(mRNA.list[[i]][mRNA.list[[i]]$pvalBH < 0.05,1:4])
}




significant.omics <-  NULL
significant.omics <- mRNA.list[[1]][mRNA.list[[1]]$pvalBH < 0.05,1:4]
significant.omics$TRAP <- "NO2"
significant.omics <- rbind(significant.omics,mRNA.list[[2]][mRNA.list[[2]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "mRNA"
significant.omics <- rbind(significant.omics,mRNA.list[[3]][mRNA.list[[3]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PM25"
significant.omics <- rbind(significant.omics,mRNA.list[[4]][mRNA.list[[4]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PCNT"
significant.omics <- rbind(significant.omics,mRNA.list[[5]][mRNA.list[[5]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "CBLK"




write.table(significant.omics, "significant_mRNAs.txt", col.names = T, row.names = FALSE)
write.csv(significant.omics, "significant_mRNAs.csv")



plots.list <- list()

for (i in 1:length(exposures.mRNA)){
  plots.list[[i]] <- ggplot(data=mRNA.list[[i]], aes(x=beta, y=-log10(pval), col = -log10(pvalBH) > -log10(0.05))) +
    geom_point(alpha=0.6, size=1) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(pvalBH) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.5) +
    geom_hline(yintercept=-log10(0.05 / nrow(mRNA.list[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot for GLS model of mRNA against", ListExpo2[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(mRNA.list[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.9) 
}

png("mRNA_volcanos.png",width = 8000, height = 2500, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[3]],plots.list[[4]],plots.list[[5]], ncol=5)
dev.off()





