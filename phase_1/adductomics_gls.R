rm(list=ls())
work_dir<-"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"

# Switch out work_dir to run on HPC or home

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"
setwd(work_dir)

#load libraries
library(mvtnorm)
library(nlme)
library(doBy)
library(parallel)


#Read imputed peaks
load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")

#Read adductomics
adductomics <- read.csv("Data/Adductomics_Oxford_Street_all_imputed.csv")

#Read Covariates and exposures 
Full_covar_add <- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
                             sep = "\t", header = T)

# Merge information from the two datasets
# Make XP_ID the rownames
rownames(Full_covar_add) <- Full_covar_add$XP_ID
rownames(Full_peaks_imputed) <- Full_peaks_imputed$XP_ID
# removing batch 56 from adductomics as it apears to be duplication
adductomics <- adductomics[adductomics$batch != 56,]
rownames(adductomics) <- adductomics$XP_ID



#log transforming adductomics
adductomics[,18:ncol(adductomics)] <- log(adductomics[,18:ncol(adductomics)])


# select people who have all information (50)
select.idx <- intersect(rownames(Full_covar_add),rownames(adductomics))

# just keep these individuals (50 individuals, 3 measurements each)
Full_covar_add <- Full_covar_add[select.idx,]
adductomics <- adductomics[select.idx,]
Full_peaks_imputed <- Full_peaks_imputed[select.idx,]

# store variables of interest
data <- NULL
data <- Full_covar_add[,c("ID","XP_ID","Group","Loc","Per","Age", 
                          "sex", "BMI", "Group", "Caffeine")]
data <- cbind(data,Full_peaks_imputed[,c("Site","Time_Point")])
data$wave <- factor(as.numeric(factor(10*as.numeric(factor(data$Site))+
                                        as.numeric(factor(data$Time_Point)))))

# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")
# Omics 
ListOMICs2 <- colnames(Full_peaks_imputed[7:5800])
ListOmics1 <- colnames(adductomics[,18:ncol(adductomics)])

# limits (repeat with FALSE, 0.05, 0.1, 0.2)
capOMIC <- TRUE
capValOMIC <- 0.05


data$Age <- as.numeric(data$Age)
data$Caffeine <- as.factor(data$Caffeine)


#Recode Group variable so that healthy becomes the reference
data$group2[data$Group == "HEALTHY"] <- "1HEALTHY"
data$group2[data$Group == "COPD"] <- "2COPD"
data$group2[data$Group == "IHD"] <- "3IHD"


##################################################
##################################################


######## FORMULA #########




NOmics <- length(ListOmics1)

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
      MyOmic <- ListOmics1[CurrOmic]
      print(paste("OMICS signal ", CurrOmic,"/",length(ListOmics1)," - ", MyOmic,sep=''))
      eval(parse(text=paste("data$Omic <- adductomics$",MyOmic,sep='')))
      
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
      
      M6 <- gls(I(Omic) ~ Expo_annual + Expo.delta + Age + sex + BMI + 
                  group2 + Caffeine,
                data=data[complete.cases(data),],
                correlation=corSymm(form=~1|ID),
                weights=varIdent(form=~1|wave),
                na.action=na.pass,
                control = list(singular.ok = TRUE))
      
      
      
      
      #Storing Results: coeffs and p-values separately
      if(MyOmic==ListOmics1[1]){
        
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
    add_resM6Pval <- cbind(ListOmics1[1:NOmics],resM6Pval)
    colnames(add_resM6Pval)[1] <- "OMIC_ID"
    add_resM6Coeff <- cbind(ListOmics1[1:NOmics],resM6Coeff)
    colnames(add_resM6Coeff)[1] <- "OMIC_ID"
    # save results 
    write.table(add_resM6Pval, paste(FileName,"_add_M6_Pvals.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    write.table(add_resM6Coeff, paste(FileName, "_add_M6_Coeff.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    
    
    
  }
}

##################################################
##################################################


##### Executing loop ####


fun.loop(ListExpo2)


##################################################
##################################################

##### Importing results #####



# read results for each exposure

add_CBLK_imp <-read.table("Results_CBLK_add_M6_Pvals.txt", header = TRUE)
add_NO2_imp <-read.table("Results_NO2_add_M6_Pvals.txt", header = TRUE)
add_PM10_imp <-read.table("Results_PM10_add_M6_Pvals.txt", header = TRUE)
add_PM25_imp <-read.table("Results_PM25_add_M6_Pvals.txt", header = TRUE)
add_PCNT_imp <-read.table("Results_PCNT_add_M6_Pvals.txt", header = TRUE)

# select significant results
add_sign_CBLK_imp <- add_CBLK_imp[(add_CBLK_imp$Expo.delta < 0.05/nrow(add_CBLK_imp)),]
add_sign_NO2_imp <- add_NO2_imp[(add_NO2_imp$Expo.delta < 0.05/nrow(add_NO2_imp)),]
add_sign_PM10_imp <- add_PM10_imp[(add_PM10_imp$Expo.delta < 0.05/nrow(add_PM10_imp)),]
add_sign_PM25_imp <- add_PM25_imp[(add_PM25_imp$Expo.delta < 0.05/nrow(add_PM25_imp)),]
add_sign_PCNT_imp <- add_PCNT_imp[(add_PCNT_imp$Expo.delta < 0.05/nrow(add_PCNT_imp)),]



##################################################
##################################################

work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/Results/phase_1/results_adductomics/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"


setwd(work_dir)

####Visualisations ####

library(ggplot2)
library(ggrepel)


# read results for each addcription exposure# 
exposures.add <- c("add_NO2_imp","add_PM10_imp","add_PM25_imp","add_PCNT_imp", "add_CBLK_imp" )


add_NO2_imp <-cbind(read.table("Results_NO2_add_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_NO2_add_M6_Coeff.txt", header = TRUE)[4])
add_PM10_imp <-cbind(read.table("Results_PM10_add_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM10_add_M6_Coeff.txt", header = TRUE)[4])
add_PM25_imp <-cbind(read.table("Results_PM25_add_M6_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM25_add_M6_Coeff.txt", header = TRUE)[4])
add_PCNT_imp <-cbind(read.table("Results_PCNT_add_M6_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_PCNT_add_M6_Coeff.txt", header = TRUE)[3])
add_CBLK_imp <-cbind(read.table("Results_CBLK_add_M6_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_CBLK_add_M6_Coeff.txt", header = TRUE)[3])

colnames(add_CBLK_imp) <- c("OMIC.ID","pval","beta")
colnames(add_NO2_imp) <- c("OMIC.ID","pval","beta")
colnames(add_PM10_imp) <- c("OMIC.ID","pval","beta")
colnames(add_PM25_imp) <- c("OMIC.ID","pval","beta")
colnames(add_PCNT_imp) <- c("OMIC.ID","pval","beta")

# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")



# Creating alist of data frames for all significant exposures
add.list <- list(add_NO2_imp,add_PM10_imp,add_PM25_imp,add_PCNT_imp,add_CBLK_imp)
names(add.list) <- exposures.add




##############################################################################################
##############################################################################################

### Creating full data frame pf pvals for manhattan plots ####

add.full.sig <- cbind(add.list[[1]][1:2], add.list[[2]][2],add.list[[3]][2], add.list[[4]][2], add.list[[5]][2])

colnames(add.full.sig) <- c("OMICid","NO2_pval","PM10_pval","PM25_pval","PCNT_pval","CBLK_pval")
write.csv(add.full.sig, "adductomics_all_pvals.csv")

##############################################################################################
##############################################################################################





for (i in 1:length(exposures.add)){
  add.list[[i]]$pvalBH <- p.adjust(add.list[[i]]$pval, method = "BH")
  add.list[[i]]$pvalFDR <- p.adjust(add.list[[i]]$pval, method = "fdr")
}


#combining data frame
significant.omics <-  NULL
significant.omics <- add.list[[1]][add.list[[1]]$pvalBH < 0.05,1:4]
significant.omics$TRAP <- "NO2"
significant.omics <- rbind(significant.omics,add.list[[2]][add.list[[2]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "miRNA"
significant.omics <- rbind(significant.omics,add.list[[3]][add.list[[3]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PM25"
significant.omics <- rbind(significant.omics,add.list[[4]][add.list[[4]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PCNT"
significant.omics <- rbind(significant.omics,add.list[[5]][add.list[[5]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "CBLK"


#write significant results to csv/txt

write.table(significant.omics, "significant_adductomics_denoised_cap05.txt", col.names = T, row.names = FALSE)
write.csv(significant.omics, "significant_adductomics_denoised_cap05.csv")


#create and save volcano plots

plots.list <- list()

for (i in 1:length(exposures.add)){
  plots.list[[i]] <- ggplot(data=add.list[[i]], aes(x=beta, y=-log10(pval), col = -log10(pvalBH) > -log10(0.05))) +
    geom_point(alpha=0.8, size=2.5) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(pvalBH) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=2.2, alpha = 0.7) +
    geom_hline(yintercept=-log10(0.05 / nrow(add.list[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot for GLS model of adductomic against", ListExpo2[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(add.list[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.9)  +
    ylim(0,-log10(0.0001))
}


png("adductomics_volcanos.png",width = 8000, height = 2500, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[3]],plots.list[[4]],plots.list[[5]], ncol=5)
dev.off()
