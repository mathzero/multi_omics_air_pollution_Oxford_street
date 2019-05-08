rm(list=ls())
setwd("~/Dropbox/Shared_HDA/Translational_Data_Science/Project6/")

#load libraries
library(mvtnorm)
library(nlme)
library(doBy)


#Read imputed peaks
load("Data/Oxford_Street_Full_peaks_imputed_metabolomics.RData")
#Read Covariates and exposures 
Full_covar_add <- 
  read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
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
capValOMIC <- 0.2


data$Age <- as.numeric(data$Age)
data$Caffeine <- as.factor(data$Caffeine)


#Recode Group variable so that healthy becomes the reference
data$group2[data$Group == "HEALTHY"] <- "1HEALTHY"
data$group2[data$Group == "COPD"] <- "2COPD"
data$group2[data$Group == "IHD"] <- "3IHD"


NOmics <- length(ListOMICs2)

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
  write.table(resM6Pval, paste(FileName,"_M6_Pvals.txt", sep = ''),
              sep = "\t", col.names = T, row.names = FALSE)
  write.table(resM6Coeff, paste(FileName, "_M6_Coeff.txt", sep = ''),
              sep = "\t", col.names = T, row.names = FALSE)
  
  
  
}
}

# parallel computing 
library(doParallel)
cl = makeCluster(detectCores())
clusterExport(cl,ls())
clusterExport(cl, c('gls', "corSymm", 'varIdent'))

system.time(parLapply(cl, ListExpo2, fun.loop))

stopCluster(cl)

# read results for each exposure
setwd('~/Desktop/Imperial/OxfordStreet/Results/nocorr/add6')
CBLK_imp <-read.table("Results_CBLK_M6_Pvalsnocorr.txt", header = TRUE)
NO2_imp <-read.table("Results_NO2_M6_Pvalsnocorr.txt", header = TRUE)
PM10_imp <-read.table("Results_PM10_M6_Pvalsnocorr.txt", header = TRUE)
PM25_imp <-read.table("Results_PM25_M6_Pvalsnocorr.txt", header = TRUE)
PCNT_imp <-read.table("Results_PCNT_M6_Pvalsnocorr.txt", header = TRUE)

# select significant results
sign_CBLK_imp <- CBLK_imp[(CBLK_imp$Expo.delta < 0.05/5794),]
sign_NO2_imp <- NO2_imp[(NO2_imp$Expo.delta < 0.05/5794),]
sign_PM10_imp <- PM10_imp[(PM10_imp$Expo.delta < 0.05/5794),]
sign_PM25_imp <- PM25_imp[(PM25_imp$Expo.delta < 0.05/5794),]
sign_PCNT_imp <- PCNT_imp[(PCNT_imp$Expo.delta < 0.05/5794),]
