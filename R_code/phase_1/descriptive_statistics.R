work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)

#Read Covariates and exposures 
Full_covar_add<- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt", sep="\t", header=T)


#plotting density graphs for exposure

#No2_r
png("NO2_levels_after1.png")
plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$NO2_r),]
plot(density(plot_subset[plot_subset$Loc == 'hp', 'NO2_r'], adjust=1), col=" dark green", main="NO2", xlab="Thousands/cm^3", ylab="Density")
lines(density(plot_subset[plot_subset$Loc == 'ox', 'NO2_r'], adjust=1), col="blue", lty=2)
legend(30,0.06, legend=c("Hyde Park", "Oxford Street"),
       col=c("dark green", "blue"), lty=1:2, cex=0.8)
dev.off()


#PM_10_r
png("PM10_levels_after1.png")
plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$PM10_r),]
plot(density(plot_subset[plot_subset$Loc == 'hp', 'PM10_r'], adjust=0.8), col=" dark green", main= "PM10", xlab="Thousands/cm^3", ylab="Density")
lines(density(plot_subset[plot_subset$Loc == 'ox', 'PM10_r'], adjust=0.8), col="blue", lty=2)
legend(60,0.04, legend=c("Hyde Park", "Oxford Street"),
       col=c("dark green", "blue"), lty=1:2, cex=0.8)
dev.off()

#PM_2.5_r
png("PM25_levels_after1.png")
plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$PM25_r),]
plot(density(plot_subset[plot_subset$Loc == 'hp', 'PM25_r'], adjust=0.8), col=" dark green", main= "PM2.5", xlab="Thousands/cm^3", ylab="Density")
lines(density(plot_subset[plot_subset$Loc == 'ox', 'PM25_r'], adjust=0.8), col="blue", lty=2)
legend(50,0.1, legend=c("Hyde Park", "Oxford Street"),
       col=c("dark green", "blue"), lty=1:2, cex=0.8)
dev.off()


# CBLK_r
png("CBLK_levels_after1.png")
plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$CBLK_r),]
plot(density(plot_subset[plot_subset$Loc == 'hp', 'CBLK_r'], adjust=1), col=" dark green", main="CBLK", xlab="Thousands/cm^3", ylab="Density")
lines(density(plot_subset[plot_subset$Loc == 'ox', 'CBLK_r'], adjust=1), col="blue", lty=2)
legend(4.5,0.3, legend=c("Hyde Park", "Oxford Street"),
       col=c("dark green", "blue"), lty=1:2, cex=0.8)
dev.off()


# PCNT_r
png("PCNT_levels_after1.png")
plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$PCNT_r),]
plot(density(plot_subset[plot_subset$Loc == 'hp', 'PCNT_r'], adjust=1), col=" dark green", main="PCNT", xlab="Thousands/cm^3", ylab="Density")
lines(density(plot_subset[plot_subset$Loc == 'ox', 'PCNT_r'], adjust=1), col="blue", lty=2)
legend(14000,0.00015, legend=c("Hyde Park", "Oxford Street"),
       col=c(" dark green", "blue"), lty=1:2, cex=0.8)
dev.off()




########################################################################
########################################################################


#plot table 1 using tableone package

#Install tableone package
if(!require(tableone)){
  install.packages("tableone")
  library(tableone)
}

library(tableone)
library(dplyr)
library(tidyverse)


#change character variables to factor (needed for categorical variables)
Full_covar_add$Loc<-as.factor(Full_covar_adds$Loc)
Full_covar_add$sex<-as.factor(Full_covar_add$sex)
Full_covar_add$season_HP<-as.factor(Full_covar_add$season_HP)
Full_covar_add$season_OS<-as.factor(Full_covar_add$season_OS)
Full_covar_add$Caffeine<-as.factor(Full_covar_add$Caffeine)

#create dataset with single observation per person
Full_covar_add2<- Full_covar_add %>% distinct(Full_covar_add$id, .keep_all = TRUE)

#Create a variable list which we want in Table 1
listVars <- c("age","Group", "BMI", "SBP", "DBP")

#create table 1
table1 <-CreateTableOne(vars = listVars, data = Full_covar_add2, strata="sex")
print(table1)

write.csv(table1, "table1.csv")

