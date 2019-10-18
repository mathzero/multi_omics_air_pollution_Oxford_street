
rm(list=ls())
### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"


#("~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/")
#"/rds/general/user/mw418/home/TDS_Project"



setwd(work_dir)

#load libraries
library(mvtnorm)
library(nlme)
library(doBy)
library(parallel)



#Read imputed peaks
load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")
#Read Covariates and exposures 
Full_covar_add <- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
                             sep = "\t", header = T)

# Merge information from the two datasets
# Make XP_ID the rownames
rownames(Full_covar_add) <- Full_covar_add$XP_ID
rownames(Full_peaks_imputed) <- Full_peaks_imputed$XP_ID

# select people who have all information (50)
select.idx <- intersect(rownames(Full_covar_add),rownames(Full_peaks_imputed))

# just keep these individuals (50 individuals, 3 measurements each)
Full_covar_add <- Full_covar_add[select.idx,]
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

# limits (repeat with FALSE, 0.05, 0.1, 0.2)
capOMIC <- TRUE
capValOMIC <- 0.05


data$Age <- as.numeric(data$Age)
data$Caffeine <- as.factor(data$Caffeine)


#Recode Group variable so that healthy becomes the reference
data$group2[data$Group == "HEALTHY"] <- "1HEALTHY"
data$group2[data$Group == "COPD"] <- "2COPD"
data$group2[data$Group == "IHD"] <- "3IHD"



###############################################################
###############################################################



NOmics <- 10 #length(ListOMICs2)

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
      MyOmic <- ListOMICs2[CurrOmic]
      print(paste("OMICS signal ", CurrOmic,"/",length(ListOMICs2)," - ", MyOmic,sep=''))
      eval(parse(text=paste("data$Omic <- Full_peaks_imputed$",MyOmic,sep='')))
      
      #standardise OMIC measurement
      data$Omic <- data$Omic/sd(data$Omic)
      
      # correct extreme values
      # select values outside quantile region and set them to the extreme 
      # repeat for upper and lower quantile
      if(capOMIC){
        data$Omic[data$Omic>quantile(data$Omic,(1-capValOMIC))]=quantile(data$Omic,(1-capValOMIC)) 
        data$Omic[data$Omic<quantile(data$Omic,capValOMIC)]=quantile(data$Omic,capValOMIC)
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
      if(MyOmic==ListOMICs2[1]){
        
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
    resM6Pval <- cbind(ListOMICs2[1:NOmics],resM6Pval)
    colnames(resM6Pval)[1] <- "OMIC_ID"
    resM6Coeff <- cbind(ListOMICs2[1:NOmics],resM6Coeff)
    colnames(resM6Coeff)[1] <- "OMIC_ID"
    # save results 
    write.table(resM6Pval, paste(FileName,"_M6_metab_Pvals.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    write.table(resM6Coeff, paste(FileName, "_M6_metab_Coeff.txt", sep = ''),
                sep = "\t", col.names = T, row.names = FALSE)
    
    
    
  }
}

###############################################################
###############################################################

### Function ####

fun.loop(ListExpo2)




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


###############################################################
###############################################################

work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/Results/phase_1/results_metabolomics/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"


setwd(work_dir)

####Visualisations ####

library(ggplot2)
library(ggrepel)


# read results for each miRNAcription exposure# 
exposures.miRNA <- c("metab_NO2_imp","metab_PM10_imp","metab_PM25_imp","metab_PCNT_imp","metab_CBLK_imp")

metab_CBLK_imp <-cbind(read.table("Results_CBLK_M6_metab_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_CBLK_M6_metab_Coeff.txt", header = TRUE)[3])
metab_NO2_imp <-cbind(read.table("Results_NO2_M6_metab_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_NO2_M6_metab_Coeff.txt", header = TRUE)[4])
metab_PM10_imp <-cbind(read.table("Results_PM10_M6_metab_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM10_M6_metab_Coeff.txt", header = TRUE)[4])
metab_PM25_imp <-cbind(read.table("Results_PM25_M6_metab_Pvals.txt", header = TRUE)[c(1,4)], read.table("Results_PM25_M6_metab_Coeff.txt", header = TRUE)[4])
metab_PCNT_imp <-cbind(read.table("Results_PCNT_M6_metab_Pvals.txt", header = TRUE)[c(1,3)], read.table("Results_PCNT_M6_metab_Coeff.txt", header = TRUE)[3])

colnames(metab_CBLK_imp) <- c("OMIC.ID","pval","beta")
colnames(metab_NO2_imp) <- c("OMIC.ID","pval","beta")
colnames(metab_PM10_imp) <- c("OMIC.ID","pval","beta")
colnames(metab_PM25_imp) <- c("OMIC.ID","pval","beta")
colnames(metab_PCNT_imp) <- c("OMIC.ID","pval","beta")


# Exposure measures
ListExpo2 <- c("NO2", "PM10", "PM25", "PCNT", "CBLK")

metab.list <- list(metab_NO2_imp,metab_PM10_imp,metab_PM25_imp,metab_PCNT_imp,metab_CBLK_imp)
names(metab.list) <- exposures.miRNA



##############################################################################################
##############################################################################################

### Creating full data frame pf pvals for manhattan plots ####

metab.full.sig <- cbind(metab.list[[1]][1:2], metab.list[[2]][2],metab.list[[3]][2], metab.list[[4]][2], metab.list[[5]][2])

colnames(metab.full.sig) <- c("OMICid","NO2_pval","PM10_pval","PM25_pval","PCNT_pval","CBLK_pval")
write.csv(metab.full.sig, "metabolites_all_pvals.csv")

##############################################################################################
##############################################################################################


for (i in 1:length(exposures.miRNA)){
  metab.list[[i]]$pvalBH <- round(p.adjust(metab.list[[i]]$pval, method = "BH"),3)
  metab.list[[i]]$pvalFDR <- round(p.adjust(metab.list[[i]]$pval, method = "fdr"),3)

  print(exposures.miRNA[i])
  print(metab.list[[i]][metab.list[[i]]$pvalBH < 0.05,1:4])
}


significant.omics <-  NULL
significant.omics <- metab.list[[1]][metab.list[[1]]$pvalBH < 0.05,1:4]
significant.omics$TRAP <- "NO2"
significant.omics <- rbind(significant.omics,metab.list[[2]][metab.list[[2]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "miRNA"
significant.omics <- rbind(significant.omics,metab.list[[3]][metab.list[[3]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PM25"
significant.omics <- rbind(significant.omics,metab.list[[4]][metab.list[[4]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "PCNT"
significant.omics <- rbind(significant.omics,metab.list[[5]][metab.list[[5]]$pvalBH < 0.05,1:4])
significant.omics$TRAP[is.na(significant.omics$TRAP)] <- "CBLK"



write.table(significant.omics, "significant_metabolites.txt", col.names = T, row.names = FALSE)
write.csv(significant.omics, "significant_metabolites.csv")


plots.list <- list()

for (i in 1:length(exposures.miRNA)){
  plots.list[[i]] <- ggplot(data=metab.list[[i]], aes(x=beta, y=-log10(pval), col = -log10(pval) > -log10(0.05 / nrow(metab.list[[i]])))) +
    geom_point(alpha=0.7, size=1) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(pval) > -log10(0.05 / nrow(metab.list[[i]])),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=2.2, alpha = 0.7) +
    geom_hline(yintercept=-log10(0.05 / nrow(metab.list[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot for GLS model of metabolomics against", ListExpo2[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(metab.list[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.6)  +
    ylim(0,-log10(0.0000000001))
}


png("metabolites_volcanos.png",width = 8000, height = 2500, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[3]],plots.list[[4]],plots.list[[5]], ncol=5)
dev.off()
