work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)


#Read Covariates and exposures 
getwd()
Full_covar_add<-read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp.txt", sep="\t", header=T)
View(Full_covar_add)


#plotting density graphs for exposure, page 29


###################################graphs###############################################################
library(ggplot2)

png(file="density_plots_exposure_combined.png",
    width=2000, height=3000,res=300)

par(mfrow = c(3, 2))

#No2_r


plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]

plot_subset <- plot_subset[!is.na(plot_subset$NO2_r),]
lim <- max(plot_subset[, 'NO2_r'])
plot(density(plot_subset[plot_subset$Loc == 'hp', 'NO2_r'], adjust=1), col=" dark green", main="NO2", xlab="Thousands/cm^3", ylab="Density", xlim=c(0,lim))
lines(density(plot_subset[plot_subset$Loc == 'ox', 'NO2_r'], adjust=1), col="blue", lty=2)
legend("topright", 
       legend = c("Hyde Park", "Oxford Street"),
       col = c("dark green","blue"),
       lty = c(1,2))


#PM_10_r

plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$PM10_r),]
lim <- max(plot_subset[, 'PM10_r'])
plot(density(plot_subset[plot_subset$Loc == 'hp', 'PM10_r'], adjust=0.8), col=" dark green", main= "PM10", xlab="Thousands/cm^3", ylab="Density", xlim=c(0,lim))
lines(density(plot_subset[plot_subset$Loc == 'ox', 'PM10_r'], adjust=0.8), col="blue", lty=2)
legend("topright", 
       legend = c("Hyde Park", "Oxford Street"),
       col = c("dark green","blue"),
       lty = c(1,2))

#PM_2.5_r


plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]

plot_subset <- plot_subset[!is.na(plot_subset$PM25_r),]
lim <- max(plot_subset[, 'PM25_r'])
plot(density(plot_subset[plot_subset$Loc == 'hp', 'PM25_r'], adjust=0.8), col=" dark green", main= "PM2.5", xlab="Thousands/cm^3", ylab="Density", xlim=c(0,lim))
lines(density(plot_subset[plot_subset$Loc == 'ox', 'PM25_r'], adjust=0.8), col="blue", lty=2)
legend("topright", 
       legend = c("Hyde Park", "Oxford Street"),
       col = c("dark green","blue"),
       lty = c(1,2))

# CBLK_r



plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$CBLK_r),]
lim <- max(plot_subset[, 'CBLK_r'])
plot(density(plot_subset[plot_subset$Loc == 'hp', 'CBLK_r'], adjust=0.8), col=" dark green", main="CBLK", xlab="Thousands/cm^3", ylab="Density", xlim=c(0,lim))
lines(density(plot_subset[plot_subset$Loc == 'ox', 'CBLK_r'], adjust=0.8), col="blue", lty=2)
legend("topright", 
       legend = c("Hyde Park", "Oxford Street"),
       col = c("dark green","blue"),
       lty = c(1,2))



# PCNT_r


plot_subset = Full_covar_add[which(Full_covar_add$Per == c('after_1')),]
plot_subset <- plot_subset[!is.na(plot_subset$PCNT_r),]
lim <- max(plot_subset[, 'PCNT_r'])
plot(density(plot_subset[plot_subset$Loc == 'hp', 'PCNT_r'], adjust=1), col=" dark green", main="PCNT", xlab="Thousands/cm^3", ylab="Density", xlim=c(0,lim))
lines(density(plot_subset[plot_subset$Loc == 'ox', 'PCNT_r'], adjust=1), col="blue", lty=2)
legend("topright", 
       legend = c("Hyde Park", "Oxford Street"),
       col = c("dark green","blue"),
       lty = c(1,2))
dev.off()



###################table one########################



#plot table 1 using table one function

library(tableone)
library(dplyr)
library(tidyverse)


#change the factors into characters (needed for categorical variables)

Full_covar_add$Loc<-as.factor(Full_covar_adds$Loc)
Full_covar_add$sex<-as.factor(Full_covar_add$sex)
Full_covar_add$season_HP<-as.factor(Full_covar_add$season_HP)
Full_covar_add$season_OS<-as.factor(Full_covar_add$season_OS)
Full_covar_add$Caffeine<-as.factor(Full_covar_add$Caffeine)
str(Full_covar_add)

#create dataset with duplicates of observations removed
Full_covar_add2<- Full_covar_add %>% distinct(Full_covar_add$id, .keep_all = TRUE)


#Create a variable list which we want in Table 1

#We need the list of ALL variables we want (in your initial script you had just included the continuous ones in listVars)
listVars <- c("age", "sex","Group", "BMI", "SBP", "DBP", "NOIS", "RHUM","temp", "NO2_r" , "PM25_r", "PM10_r", "CBLK_r","PCNT_r")
#Factors
catVars <- c("sex", "Group")


#Note that since you made sure all your factors are actually factors beforehand,
# factorVars argument isn't necessary
table1 <-CreateTableOne( vars = listVars, data = Full_covar_add2, strata="sex")
print(table1)

#how to further stratify for males and females in group? that is to be done
table2 <-CreateTableOne( vars = listVars, data = Full_covar_add2, strata ="Group", "sex")
print(table2)





########################pheatmaps##################################

#pheatmap amongst exposures (before walk- create dataframe with just before walk)
library(pheatmap)

dfbefore<-Full_covar_add[which(Full_covar_add$Per=='before'),]

#for dfbefore select exposures in dataframe
dfbefore1<-data.frame(dfbefore$NO2_r, dfbefore$PM10_r, dfbefore$PM25_r, dfbefore$PCNT_r, dfbefore$CBLK_r)
colnames(dfbefore1)
colnames(dfbefore1)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(before walk)
mycor=cor(dfbefore1, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresbefore.png",
    width=600, height=350)
pheatmap(mycor, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures Before")
dev.off()



#pheat map amongst exposures (after walk- create dataframe with just after 1)
dfafter1<- Full_covar_add[which(Full_covar_add$Per=='after_1'),]


#for dfbefore select exposures in dataframe after 1
dfafter1x<-data.frame(dfafter1$NO2_r, dfafter1$PM10_r, dfafter1$PM25_r, dfafter1$PCNT_r, dfafter1$CBLK_r)
colnames(dfafter1x)
colnames(dfafter1x)= c("NO2","PM10","PM25","PCNT","CBLK")

#pheat map of exposures_r(after walk 1)
mycor1=cor(dfafter1x, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresafter1.png",
    width=600, height=350)
pheatmap(mycor1, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures after 1")
dev.off()



#pheat map amongst exposures (after walk- create dataframe with just after 2)
dfafter2<- Full_covar_add[which(Full_covar_add$Per=='after_2'),]

#for dfbefore select exposures in dataframe after 2
dfafter2x<-data.frame(dfafter2$NO2_r, dfafter2$PM10_r, dfafter2$PM25_r, dfafter2$PCNT_r, dfafter2$CBLK_r)
colnames(dfafter2x)
colnames(dfafter2x)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(after walk 2)
mycor2=cor(dfafter2x, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresafter2.png",
    width=600, height=350)
pheatmap(mycor2, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures after 2")
dev.off()



#pheat map amongst exposures (delta -not really useful just have it as a test)

#create a new dataframe with delta exposures as variables- use the preprogrammed deltas- I am just going to do everything- annual

Full_covar_add3<-Full_covar_add
Full_covar_add3$NO2delta<-Full_covar_add3$NO2_r- Full_covar_add3$NO2_annual
Full_covar_add3$PM10delta<-Full_covar_add$PM10_r-Full_covar_add3$PM10_annual
Full_covar_add3$PM25delta<-Full_covar_add3$PM25_r-Full_covar_add3$PM25_annual
Full_covar_add3$PCNTdelta<-Full_covar_add3$PCNT_r-Full_covar_add3$PCNT_annual
Full_covar_add3$CBLKdelta<-Full_covar_add3$CBLK_r-Full_covar_add3$CBLK_annual

dfdeltaexpo<-data.frame(Full_covar_add3$NO2delta, Full_covar_add3$PM10delta, Full_covar_add3$PM25delta, Full_covar_add3$PCNTdelta, Full_covar_add3$CBLKdelta)


#pheat map dfdeltaexpo- exposures delta (no stratification for before or after)
mycor3=cor(dfdeltaexpo, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresdeltaannaul.png",
    width=600, height=350)
pheatmap(mycor3, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100))
dev.off()


#pheat map after2-after1
# try to merge the two dataframes
dfmerged1and2<-merge(x = dfafter1, y = dfafter2, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures
dfmerged1and2$NO2delta<-dfmerged1and2$NO2_r.y- dfmerged1and2$NO2_r.x
dfmerged1and2$PM10delta<-dfmerged1and2$PM10_r.y-dfmerged1and2$PM10_r.x
dfmerged1and2$PM25delta<-dfmerged1and2$PM25_r.y-dfmerged1and2$PM25_r.x
dfmerged1and2$PCNTdelta<-dfmerged1and2$PCNT_r.y-dfmerged1and2$PCNT_r.x
dfmerged1and2$CBLKdelta<-dfmerged1and2$CBLK_r.y-dfmerged1and2$CBLK_r.x

#create dataframe for just deltas (after2-after1)

dfdeltaexpoafter1and2<-data.frame(dfmerged1and2$NO2delta, dfmerged1and2$PM10delta, dfmerged1and2$PM25delta, dfmerged1and2$PCNTdelta, dfmerged1and2$CBLKdelta)

colnames(dfdeltaexpoafter1and2)
colnames(dfdeltaexpoafter1and2)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 2- after 1
mycor4=cor(dfdeltaexpoafter1and2, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/delta exposures (after 2-1).png",
    width=600, height=350)
pheatmap(mycor4, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after2- after1)")

dev.off()




##labeling ways : just to rememberpheatmap(mat, labels_row=paste0("foo", 1:10), labels_col=paste0("bar", 1:10))
#seperating location, and time and delta each time point, and each location
#igraph ::, intersection





#pheat map after2-before

# try to merge the two dataframes
dfmerged2andb<-merge(x = dfbefore, y = dfafter2, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures
dfmerged2andb$NO2delta<-dfmerged2andb$NO2_r.y- dfmerged2andb$NO2_r.x
dfmerged2andb$PM10delta<-dfmerged2andb$PM10_r.y-dfmerged2andb$PM10_r.x
dfmerged2andb$PM25delta<-dfmerged2andb$PM25_r.y-dfmerged2andb$PM25_r.x
dfmerged2andb$PCNTdelta<-dfmerged2andb$PCNT_r.y-dfmerged2andb$PCNT_r.x
dfmerged2andb$CBLKdelta<-dfmerged2andb$CBLK_r.y-dfmerged2andb$CBLK_r.x

#create dataframe for just deltas (after2-before)

dfdeltaexpoafter2andb<-data.frame(dfmerged2andb$NO2delta, dfmerged2andb$PM10delta, dfmerged2andb$PM25delta, dfmerged2andb$PCNTdelta, dfmerged2and$CBLKdelta)
colnames(dfdeltaexpoafter2andb)
colnames(dfdeltaexpoafter2andb)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 2- before
mycor5=cor(dfdeltaexpoafter2andb, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/deltaexposuresafter2-before.png",
    width=600, height=350)
pheatmap(mycor5, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after2- before)")
dev.off()
###

#pheat map after1-before

# try to merge the two dataframes
dfmerged1andb<-merge(x = dfbefore, y = dfafter1, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures
dfmerged1andb$NO2delta<-dfmerged1andb$NO2_r.y- dfmerged1andb$NO2_r.x
dfmerged1andb$PM10delta<-dfmerged1andb$PM10_r.y-dfmerged1andb$PM10_r.x
dfmerged1andb$PM25delta<-dfmerged1andb$PM25_r.y-dfmerged1andb$PM25_r.x
dfmerged1andb$PCNTdelta<-dfmerged1andb$PCNT_r.y-dfmerged1andb$PCNT_r.x
dfmerged1andb$CBLKdelta<-dfmerged1andb$CBLK_r.y-dfmerged1andb$CBLK_r.x

#create dataframe for just deltas (after1-before)

dfdeltaexpoafter1andb<-data.frame(dfmerged1andb$NO2delta, dfmerged1andb$PM10delta, dfmerged1andb$PM25delta, dfmerged1andb$PCNTdelta, dfmerged1andb$CBLKdelta)
colnames(dfdeltaexpoafter1andb)
colnames(dfdeltaexpoafter1andb)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 1- before
mycor6=cor(dfdeltaexpoafter1andb, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/deltaexpo-after1-before.png",
    width=600, height=350)
pheatmap(mycor6, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after1- before)")

dev.off()
#########################################################

##pheat maps seperated by location

#pheatmap amongst exposures (before walk- create dataframe with just before walk)
library(pheatmap)
dfbeforehp<-Full_covar_add[which(Full_covar_add$Per=='before' & Full_covar_add$Loc=='hp'),]


#for dfbefore select exposures in dataframe
dfbefore1hp<-data.frame(dfbeforehp$NO2_r, dfbeforehp$PM10_r, dfbeforehp$PM25_r, dfbeforehp$PCNT_r, dfbeforehp$CBLK_r)
colnames(dfbefore1hp)
colnames(dfbefore1hp)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(before walk hp)
mycorhp=cor(dfbefore1hp, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresbeforeHP.png",
    width=600, height=350)
pheatmap(mycorhp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures Before Walk in Hyde Park")
dev.off()

#pheat map amongst exposures (after 1 walk- create dataframe with just after 1 hp)
dfafter1hp<- Full_covar_add[which(Full_covar_add$Per=='after_1' & Full_covar_add$Loc=='hp'),]


#for dfafter select exposures in dataframe after 1 hp
dfafter1xhp<-data.frame(dfafter1hp$NO2_r, dfafter1hp$PM10_r, dfafter1hp$PM25_r, dfafter1hp$PCNT_r, dfafter1hp$CBLK_r)
colnames(dfafter1xhp)
colnames(dfafter1xhp)= c("NO2","PM10","PM25","PCNT","CBLK")

#pheat map of exposures_r(after walk 1 hp)
mycor1hp=cor(dfafter1xhp, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresafter1hp.png",
    width=600, height=350)
pheatmap(mycor1hp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures 'after 1' in Hyde Park")
dev.off()

#pheat map of exposure_r (after walk 2 hp)

#pheat map amongst exposures (after walk- create dataframe with just after 2 hp)
dfafter2hp<- Full_covar_add[which(Full_covar_add$Per=='after_2'& Full_covar_add$Loc=='hp'),]


#for dfbefore select exposures in dataframe after 2-hyde park
dfafter2xhp<-data.frame(dfafter2hp$NO2_r, dfafter2hp$PM10_r, dfafter2hp$PM25_r, dfafter2hp$PCNT_r, dfafter2hp$CBLK_r)
colnames(dfafter2xhp)
colnames(dfafter2xhp)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(after walk 2 hyde park)
mycor2hp=cor(dfafter2xhp, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresafter2hp.png",
    width=600, height=350)
pheatmap(mycor2hp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures 'after 2' in Hyde Park")
dev.off()

##########pheatmaps for delta exposure by location########################


#pheat map after2-after1 hyde park

# try to merge the two dataframes
dfmerged1and2hp<-merge(x = dfafter1hp, y = dfafter2hp, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures
dfmerged1and2hp$NO2delta<-dfmerged1and2hp$NO2_r.y- dfmerged1and2hp$NO2_r.x
dfmerged1and2hp$PM10delta<-dfmerged1and2hp$PM10_r.y-dfmerged1and2hp$PM10_r.x
dfmerged1and2hp$PM25delta<-dfmerged1and2hp$PM25_r.y-dfmerged1and2hp$PM25_r.x
dfmerged1and2hp$PCNTdelta<-dfmerged1and2hp$PCNT_r.y-dfmerged1and2hp$PCNT_r.x
dfmerged1and2hp$CBLKdelta<-dfmerged1and2hp$CBLK_r.y-dfmerged1and2hp$CBLK_r.x

#create dataframe for just deltas (after2-after1)

dfdeltaexpoafter1and2hp<-data.frame(dfmerged1and2hp$NO2delta, dfmerged1and2hp$PM10delta, dfmerged1and2hp$PM25delta, dfmerged1and2hp$PCNTdelta, dfmerged1and2hp$CBLKdelta)

colnames(dfdeltaexpoafter1and2hp)
colnames(dfdeltaexpoafter1and2hp)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 2- after 1
mycor4hp=cor(dfdeltaexpoafter1and2hp, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/after2-1hp.png",
    width=600, height=350)
pheatmap(mycor4hp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after2- after1), hyde park")

dev.off()

##labeling ways : just to rememberpheatmap(mat, labels_row=paste0("foo", 1:10), labels_col=paste0("bar", 1:10))
#seperating location, and time and delta each time point, and each location
#igraph ::, intersection





#pheat map after2-before (hyde park)

# try to merge the two dataframes
dfmerged2andbhp<-merge(x = dfbeforehp, y = dfafter2hp, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures- hyde park
dfmerged2andbhp$NO2delta<-dfmerged2andbhp$NO2_r.y- dfmerged2andbhp$NO2_r.x
dfmerged2andbhp$PM10delta<-dfmerged2andbhp$PM10_r.y-dfmerged2andbhp$PM10_r.x
dfmerged2andbhp$PM25delta<-dfmerged2andbhp$PM25_r.y-dfmerged2andbhp$PM25_r.x
dfmerged2andbhp$PCNTdelta<-dfmerged2andbhp$PCNT_r.y-dfmerged2andbhp$PCNT_r.x
dfmerged2andbhp$CBLKdelta<-dfmerged2andbhp$CBLK_r.y-dfmerged2andbhp$CBLK_r.x

#create dataframe for just deltas (after2-before- hyde park)

dfdeltaexpoafter2andbhp<-data.frame(dfmerged2andbhp$NO2delta, dfmerged2andbhp$PM10delta, dfmerged2andbhp$PM25delta, dfmerged2andbhp$PCNTdelta, dfmerged2andbhp$CBLKdelta)

colnames(dfdeltaexpoafter2andbhp)
colnames(dfdeltaexpoafter2andbhp)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 2- before- hyde park
mycor5hp=cor(dfdeltaexpoafter2andbhp, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/after2-beforehp.png",
    width=600, height=350)
pheatmap(mycor5hp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after2- before), hyde park")
dev.off()
###

#pheat map after1-before- hyde park

# try to merge the two dataframes
dfmerged1andbhp<-merge(x = dfbeforehp, y = dfafter1hp, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures
dfmerged1andbhp$NO2delta<-dfmerged1andbhp$NO2_r.y- dfmerged1andbhp$NO2_r.x
dfmerged1andbhp$PM10delta<-dfmerged1andbhp$PM10_r.y-dfmerged1andbhp$PM10_r.x
dfmerged1andbhp$PM25delta<-dfmerged1andbhp$PM25_r.y-dfmerged1andbhp$PM25_r.x
dfmerged1andbhp$PCNTdelta<-dfmerged1andbhp$PCNT_r.y-dfmerged1andbhp$PCNT_r.x
dfmerged1andbhp$CBLKdelta<-dfmerged1andb$CBLK_r.y-dfmerged1andbhp$CBLK_r.x

#create dataframe for just deltas (after1-before), hyde park

dfdeltaexpoafter1andbhp<-data.frame(dfmerged1andbhp$NO2delta, dfmerged1andbhp$PM10delta, dfmerged1andbhp$PM25delta, dfmerged1andbhp$PCNTdelta, dfmerged1andbhp$CBLKdelta)

colnames(dfdeltaexpoafter1andbhp)
colnames(dfdeltaexpoafter1andbhp)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 1- before, hyde park
mycor6hp=cor(dfdeltaexpoafter1andb, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/after1-beforehp.png",
    width=600, height=350)
pheatmap(mycor6hp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after1- before), hyde park")
dev.off()
###############oxford street location#######
##################################################################################


##pheat maps seperated by location

#pheatmap amongst exposures (before walk- create dataframe with just before walk)
library(pheatmap)
dfbeforeox<-Full_covar_add[which(Full_covar_add$Per=='before' & Full_covar_add$Loc=='ox'),]


#for dfbefore select exposures in dataframe
dfbefore1ox<-data.frame(dfbeforeox$NO2_r, dfbeforeox$PM10_r, dfbeforeox$PM25_r, dfbeforeox$PCNT_r, dfbeforeox$CBLK_r)
colnames(dfbefore1ox)
colnames(dfbefore1ox)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(before walk ox)
mycorox=cor(dfbefore1ox, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresbeforeox.png",
    width=600, height=350)
pheatmap(mycorox, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures Before Walk in Oxford Street")
dev.off()

#pheat map amongst exposures (after 1 walk- create dataframe with just after 1 ox)
dfafter1ox<- Full_covar_add[which(Full_covar_add$Per=='after_1' & Full_covar_add$Loc=='ox'),]


#for dfafter select exposures in dataframe after 1 ox
dfafter1xox<-data.frame(dfafter1ox$NO2_r, dfafter1ox$PM10_r, dfafter1ox$PM25_r, dfafter1ox$PCNT_r, dfafter1ox$CBLK_r)
colnames(dfafter1xox)
colnames(dfafter1xox)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(after walk 1 ox)
mycor1hp=cor(dfafter2xhp, use="complete.obs")
mycor1ox=cor(dfafter1xox, use="complete.obs")

png(file="after1walkox.png",
    width=2000, height=2000, res = 300)
pheatmap(mycor1ox, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposure correlation during walk in Oxford Street")
dev.off()

png(file="after1walkhp.png",
    width=2000, height=2000, res = 300)
pheatmap(mycor1hp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposure correlation during walk in Hyde Park")
dev.off()
#pheat map of exposure_r (after walk 2 ox)

#pheat map amongst exposures (after walk- create dataframe with just after 2 ox)
dfafter2ox<- Full_covar_add[which(Full_covar_add$Per=='after_2'& Full_covar_add$Loc=='ox'),]


#for dfbefore select exposures in dataframe after 2-ox
dfafter2xox<-data.frame(dfafter2ox$NO2_r, dfafter2ox$PM10_r, dfafter2ox$PM25_r, dfafter2ox$PCNT_r, dfafter2ox$CBLK_r)
colnames(dfafter2xox)
colnames(dfafter2xox)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r(after walk 2 ox)
mycor2ox=cor(dfafter2xox, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/after2ox.png",
    width=600, height=350)
pheatmap(mycor2ox, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures 'after 2' in Oxford Street")
dev.off()

##########pheatmaps for delta exposure by location########################


#pheat map after2-after1 oxford st

# try to merge the two dataframes
dfmerged1and2ox<-merge(x = dfafter1ox, y = dfafter2ox, by = "ID", all = TRUE)

#create new variables in dfmerged1and2 for delta exposures
dfmerged1and2ox$NO2delta<-dfmerged1and2ox$NO2_r.y- dfmerged1and2ox$NO2_r.x
dfmerged1and2ox$PM10delta<-dfmerged1and2ox$PM10_r.y-dfmerged1and2ox$PM10_r.x
dfmerged1and2ox$PM25delta<-dfmerged1and2ox$PM25_r.y-dfmerged1and2ox$PM25_r.x
dfmerged1and2ox$PCNTdelta<-dfmerged1and2ox$PCNT_r.y-dfmerged1and2ox$PCNT_r.x
dfmerged1and2ox$CBLKdelta<-dfmerged1and2hp$CBLK_r.y-dfmerged1and2ox$CBLK_r.x

#create dataframe for just deltas (after2-after1)

dfdeltaexpoafter1and2ox<-data.frame(dfmerged1and2ox$NO2delta, dfmerged1and2ox$PM10delta, dfmerged1and2ox$PM25delta, dfmerged1and2ox$PCNTdelta, dfmerged1and2ox$CBLKdelta)

colnames(dfdeltaexpoafter1and2ox)
colnames(dfdeltaexpoafter1and2ox)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 2- after 1
mycor4ox=cor(dfdeltaexpoafter1and2ox, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/delta2-1ox.png",
    width=600, height=350)
pheatmap(mycor4ox, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after2- after1), oxford street")

dev.off()

##labeling ways : just to rememberpheatmap(mat, labels_row=paste0("foo", 1:10), labels_col=paste0("bar", 1:10))
#seperating location, and time and delta each time point, and each location
#igraph ::, intersection





#pheat map after2-before (oxford street)

# try to merge the two dataframes
dfmerged2andbox<-merge(x = dfbeforeox, y = dfafter2ox, by = "ID", all = TRUE)


#create new variables in dfmerged1and2 for delta exposures- ox
dfmerged2andbox$NO2delta<-dfmerged2andbox$NO2_r.y- dfmerged2andbox$NO2_r.x
dfmerged2andbox$PM10delta<-dfmerged2andbox$PM10_r.y-dfmerged2andbox$PM10_r.x
dfmerged2andbox$PM25delta<-dfmerged2andbox$PM25_r.y-dfmerged2andbox$PM25_r.x
dfmerged2andbox$PCNTdelta<-dfmerged2andbox$PCNT_r.y-dfmerged2andbox$PCNT_r.x
dfmerged2andbox$CBLKdelta<-dfmerged2andbox$CBLK_r.y-dfmerged2andbox$CBLK_r.x

#create dataframe for just deltas (after2-before- ox)

dfdeltaexpoafter2andbox<-data.frame(dfmerged2andbox$NO2delta, dfmerged2andbox$PM10delta, dfmerged2andbox$PM25delta, dfmerged2andbox$PCNTdelta, dfmerged2andbox$CBLKdelta)

colnames(dfdeltaexpoafter2andbox)
colnames(dfdeltaexpoafter2andbox)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 2- before-ox
mycor5ox=cor(dfdeltaexpoafter2andbox, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/after2-box.png",
    width=600, height=350)
pheatmap(mycor5ox, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after2- before), oxford street")
dev.off()
###

#pheat map after1-before- ox

# try to merge the two dataframes
dfmerged1andbox<-merge(x = dfbeforeox, y = dfafter1ox, by = "ID", all = TRUE)


#create new variables in dfmerged1and b for delta exposures
dfmerged1andbox$NO2delta<-dfmerged1andbox$NO2_r.y- dfmerged1andbox$NO2_r.x
dfmerged1andbox$PM10delta<-dfmerged1andbox$PM10_r.y-dfmerged1andbox$PM10_r.x
dfmerged1andbox$PM25delta<-dfmerged1andbox$PM25_r.y-dfmerged1andbox$PM25_r.x
dfmerged1andbox$PCNTdelta<-dfmerged1andbox$PCNT_r.y-dfmerged1andbox$PCNT_r.x
dfmerged1andbox$CBLKdelta<-dfmerged1andbox$CBLK_r.y-dfmerged1andbox$CBLK_r.x

#create dataframe for just deltas (after1-before), ox

dfdeltaexpoafter1andbox<-data.frame(dfmerged1andbox$NO2delta, dfmerged1andbox$PM10delta, dfmerged1andbox$PM25delta, dfmerged1andbox$PCNTdelta, dfmerged1andbox$CBLKdelta)

colnames(dfdeltaexpoafter1andbox)
colnames(dfdeltaexpoafter1andbox)= c("NO2","PM10","PM25","PCNT","CBLK")


#create pheatmap for after 1- before, ox
mycor6ox=cor(dfdeltaexpoafter1andbox, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/after1-box.png",
    width=600, height=350)
pheatmap(mycor6ox, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="delta of exposures (after1- before), oxford street")
dev.off()

#########################################################################################################################################################


##########pheat maps just location#################

# 
# Oxford street
library(pheatmap)
dfox<-Full_covar_add[which( Full_covar_add$Loc=='ox'),]


#for dfbefore select exposures in dataframe
dfox<-data.frame(dfox$NO2_r, dfox$PM10_r, dfox$PM25_r, dfox$PCNT_r, dfox$CBLK_r)
colnames(dfox)
colnames(dfox)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r ox
mycoroxs=cor(dfox, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresoxfordstreet.png",
    width=600, height=350)
pheatmap(mycoroxs, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures  Oxford Street")
dev.off()

##hyde park
library(pheatmap)
dfhp<-Full_covar_add[which( Full_covar_add$Loc=='hp'),]


#for dfbefore select exposures in dataframe
dfhp<-data.frame(dfhp$NO2_r, dfhp$PM10_r, dfhp$PM25_r, dfhp$PCNT_r, dfhp$CBLK_r)
colnames(dfhp)
colnames(dfhp)= c("NO2","PM10","PM25","PCNT","CBLK")


#pheat map of exposures_r ox
mycorhp=cor(dfhp, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposureshp.png",
    width=600, height=350)
pheatmap(mycorhp, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main="Exposures Hyde Park")

dev.off()

# try to merge the two dataframes
dfdelta<-dfox-dfhp

#create pheatmap for after 2- after 1
mycordelta=cor(dfdelta, use="complete.obs")
png(file="~/Documents/R_materials_and_datasets/exposuresdeltaoxhp.png",
    width=600, height=350)
pheatmap(mycordelta, cluster_cols=FALSE, cluster_rows=FALSE, breaks=seq(-1,1, length.out=100), main=" Correlation of delta of exposures/regardless of time")

dev.off()