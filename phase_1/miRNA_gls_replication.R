rm(list=ls())
work_dir<-"/rds/general/user/mw418/home/TDS_Project"

# Switch out work_dir to run on HPC or home

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"


setwd(work_dir)
#load libraries
library(mvtnorm)
library(lme4)
library(nlme)
library(doBy)
library(dplyr)
library(parallel)


#Read imputed peaks
load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")

###### Read in miRNAomics ####

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

miRNA$XP_ID <- paste(miRNA$ID,numbers, sep = "")

miRNA <- miRNA[,c(373:375,1:372)] 


##############


#Read Covariates and exposures 
Full_covar_add <- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
                             sep = "\t", header = T)

# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")


# testing running the code with t0 as annual #####



# ###### New annual figures ####
# Full_covar_add$CBLK_annualxx <-  NA
# Full_covar_add$NO2_annualxx <-  NA
# Full_covar_add$PM10_annualxx <-  NA
# Full_covar_add$PM25_annualxx <-  NA
# Full_covar_add$PCNT_annualxx <-  NA
# 
# 
# #######   For loop to add t0 values as new 'annualxx' values #####
# 
# for (expo in ListExpo2){
#   for (i in 1:length(Full_covar_add$ID)){
#     print(Full_covar_add$XP_ID[i])
#     ID <- as.numeric(gsub('.{2}$', '', Full_covar_add$XP_ID[i]))
#     if (Full_covar_add$Loc[i] == "hp") {
#       if (length(eval(parse(text=paste0("Full_covar_add$", expo ,"_o",sep='')))[Full_covar_add$ID == ID & Full_covar_add$Per == "before" & Full_covar_add$Loc == "hp"]) == 0){
#         print(paste(i,"Oh No HP!"))
#         next
#       }
#       else{
#         eval(parse(text=paste0("Full_covar_add$", expo ,"_annualxx[i] <- eval(parse(text=paste0(\"Full_covar_add$\", expo ,\"_o\",sep='')))[Full_covar_add$ID == ID & 
#                                                                    Full_covar_add$Per == \"before\" & 
#                                                                    Full_covar_add$Loc == \"hp\"]")))
#       }
#     }
#     else{
#       if (length(eval(parse(text=paste0("Full_covar_add$", expo ,"_o",sep='')))[Full_covar_add$ID == ID & Full_covar_add$Per == "before" & Full_covar_add$Loc == "ox"]) == 0){
#         print(paste(i,"Oh No Ox!"))
#         next
#       }
#       else{
#         eval(parse(text=paste0("Full_covar_add$", expo ,"_annualxx[i] <- eval(parse(text=paste0(\"Full_covar_add$\", expo ,\"_o\",sep='')))[Full_covar_add$ID == ID & 
#                                                                    Full_covar_add$Per == \"before\" & 
#                                                                    Full_covar_add$Loc == \"ox\"] "
#         )))
#       }
#     }
#   }
#   
# }
# 

###############################################################
###############################################################


#########


# Merge information from the two datasets
# Make XP_ID the rownames
rownames(Full_covar_add) <- Full_covar_add$XP_ID
rownames(Full_peaks_imputed) <- Full_peaks_imputed$XP_ID


rownames(miRNA) <- miRNA$XP_ID




# select people who have all information (50)
select.idx <- intersect(rownames(Full_covar_add),rownames(miRNA))

# just keep these individuals (50 individuals, 3 measurements each)
Full_covar_add <- Full_covar_add[select.idx,]
miRNA <- miRNA[select.idx,]
Full_peaks_imputed <- Full_peaks_imputed[select.idx,]

# store variables of interest
data <- NULL
data <- Full_covar_add[,c("ID","XP_ID","Group","Loc","Per","Age", 
                          "sex", "BMI", "Group", "Caffeine")]
data <- cbind(data,Full_peaks_imputed[,c("Site","Time_Point")])
data <- cbind(data,miRNA[,c("Batch","Array")])
data$wave <- factor(as.numeric(factor(10*as.numeric(factor(data$Site))+
                                        as.numeric(factor(data$Time_Point)))))


#removing bizarre miRNA with same value for every observation###

miRNA <-  subset(miRNA, select=-hsa.miR.451a)


# Omics 
ListOMICs2 <- colnames(Full_peaks_imputed[7:5800])
ListOmics3 <- colnames(miRNA[,10:ncol(miRNA)])

# limits (repeat with FALSE, 0.05, 0.1, 0.2)
capOMIC <- TRUE
capValOMIC <- 0.2


data$Age <- as.numeric(data$Age)
data$Caffeine <- as.factor(data$Caffeine)


#Recode Group variable so that healthy becomes the reference
data$group2[data$Group == "HEALTHY"] <- "1HEALTHY"
data$group2[data$Group == "COPD"] <- "2COPD"
data$group2[data$Group == "IHD"] <- "3IHD"



###############################################################
###############################################################





######## FORMULA #########


NOmics <- length(ListOmics3)

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
      MyOmic <- ListOmics3[CurrOmic]
      print(paste("OMICS signal ", CurrOmic,"/",length(ListOmics3)," - ", MyOmic,sep=''))
      eval(parse(text=paste("data$Omic <- miRNA$",MyOmic,sep='')))
      
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
      
      preds <- lmer(Omic ~(1|Array) + (1|Batch), data = data, REML = TRUE)
      data$Omic <- resid(preds)
      
      M6 <- gls(I(Omic) ~ Expo_annual + Expo.delta + Age + sex + BMI + 
                  group2,
                data=data[complete.cases(data),],
                correlation=corSymm(form=~1|ID),
                weights=varIdent(form=~1|wave),
                na.action=na.pass,
                control = list(singular.ok = TRUE))
      
      
      
      
      #Storing Results: coeffs and p-values separately
      if(MyOmic==ListOmics3[1]){
        
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
    miRNA_resM6Pval <- cbind(ListOmics3[1:NOmics],resM6Pval)
    colnames(miRNA_resM6Pval)[1] <- "OMIC_ID"
    miRNA_resM6Coeff <- cbind(ListOmics3[1:NOmics],resM6Coeff)
    colnames(miRNA_resM6Coeff)[1] <- "OMIC_ID"
    # save results 
    write.table(miRNA_resM6Pval, paste(FileName,"_miRNA_M6_Pvals.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    write.table(miRNA_resM6Coeff, paste(FileName, "_miRNA_M6_Coeff.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    
    
    
  }
}




###############################################################
###############################################################

########## Executing loop ########

###### Single core (uncomment to execute - shift control C)  ######




fun.loop(ListExpo2)




######  Parrallelised (uncomment to execute)  ######
# 
# 
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
# 


###############################################################
###############################################################



####Visualisations ####

library(ggplot2)
library(ggrepel)

# read results for each miRNAcription exposure# 
exposures.miRNA <- c("miRNA_NO2_imp","miRNA_PM10_imp","miRNA_PM25_imp","miRNA_PCNT_imp","miRNA_CBLK_imp")

miRNA_NO2_imp <-cbind(read.table("Results_NO2_miRNA_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_NO2_miRNA_M6_Coeff.txt", header = TRUE)[4])
miRNA_PM10_imp <-cbind(read.table("Results_PM10_miRNA_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM10_miRNA_M6_Coeff.txt", header = TRUE)[4])
miRNA_PM25_imp <-cbind(read.table("Results_PM25_miRNA_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM25_miRNA_M6_Coeff.txt", header = TRUE)[4])
miRNA_PCNT_imp <-cbind(read.table("Results_PCNT_miRNA_M6_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_PCNT_miRNA_M6_Coeff.txt", header = TRUE)[3])
miRNA_CBLK_imp <-cbind(read.table("Results_CBLK_miRNA_M6_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_CBLK_miRNA_M6_Coeff.txt", header = TRUE)[3])


colnames(miRNA_NO2_imp) <- c("OMIC.ID","pval","beta")
colnames(miRNA_PM10_imp) <- c("OMIC.ID","pval","beta")
colnames(miRNA_PM25_imp) <- c("OMIC.ID","pval","beta")
colnames(miRNA_PCNT_imp) <- c("OMIC.ID","pval","beta")
colnames(miRNA_CBLK_imp) <- c("OMIC.ID","pval","beta")

# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")

miRNA.list <- list(miRNA_NO2_imp,miRNA_PM10_imp,miRNA_PM25_imp,miRNA_PCNT_imp,miRNA_CBLK_imp)
names(miRNA.list) <- exposures.miRNA

for (i in 1:length(exposures.miRNA)){
  miRNA.list[[i]]$pvalBH <- p.adjust(miRNA.list[[i]]$pval, method = "BH")
  miRNA.list[[i]]$pvalFDR <- p.adjust(miRNA.list[[i]]$pval, method = "fdr")
}




significant.omics <-  NULL
significant.omics <- miRNA.list[[1]][miRNA.list[[1]]$pvalBH < 0.05,1:4]
significant.omics$TRAP <- "NO2"
significant.omics <- rbind(significant.omics,miRNA.list[[2]][miRNA.list[[2]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "miRNA"
significant.omics <- rbind(significant.omics,miRNA.list[[3]][miRNA.list[[3]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PM25"
significant.omics <- rbind(significant.omics,miRNA.list[[4]][miRNA.list[[4]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PCNT"
significant.omics <- rbind(significant.omics,miRNA.list[[5]][miRNA.list[[5]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "CBLK"


write.table(significant.omics, "significant_miRNA_denoised_replication.txt", col.names = T, row.names = FALSE)
write.csv(significant.omics, "significant_miRNA_denoised_replication.csv")


for (i in 1:length(exposures.miRNA)){
  g <- ggplot(data=miRNA.list[[i]], aes(x=beta, y=-log10(pval), col = -log10(pvalBH) > -log10(0.05))) +
    geom_point(alpha=0.8, size=1) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(pvalBH) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=2.2, alpha = 0.7) +
    geom_hline(yintercept=-log10(0.05 / nrow(miRNA.list[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot for GLS model of miRNA against", ListExpo2[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(miRNA.list[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.6) 
  print(g)
  ggsave(filename = paste0(exposures.miRNA[i],"_miRNA_volcano.png"), plot = last_plot(), device = "png", 
         dpi = 300)
}



