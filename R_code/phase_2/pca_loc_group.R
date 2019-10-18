#### Machine learning for health status ####
library(ggbiplot)
library(ggfortify)
library(lfda)

#### Metabolites#####
Full_covar_add$Per <- factor(Full_covar_add$Per, levels = c("before", "after_1", "after_2"))

times <- levels(metabolites$Time_Point)

group.pca.plots <- list()
loc.pca.plots <- list()

for (i in times){
  omic.time <-  metabolites[metabolites$Time_Point == i,]
  omic.pcr <- prcomp(omic.time[,7:(ncol(omic.time)-3)], center = TRUE,scale. = TRUE)
  group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group") + ggtitle(i)
  loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
}


png("metabolites_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
dev.off()

png("metabolites_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
dev.off()


### mRNA ###

mRNA <- inner_join(Full_covar_add, mRNA, by = "XP_ID")

times <- levels(mRNA$Time_Point)

group.pca.plots <- list()
loc.pca.plots <- list()


for (i in times){
  omic.time <-  mRNA[mRNA$Time_Point == i,]
  omic.pcr <- prcomp(omic.time[,59:(ncol(omic.time))], center = TRUE,scale. = TRUE)
  group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group") + ggtitle(i)
  loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
}


png("mRNA_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
dev.off()

png("mRNA_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
dev.off()


#### miRNA #####

miRNA <- inner_join(Full_covar_add, miRNA, by = "XP_ID")

times <- levels(miRNA$Time_Point)

miRNA <- miRNA[-as.numeric(which(apply(miRNA, 2, var) == 0))] #removing columns with zero variance

miRNA[1:3,50:60]

group.pca.plots <- list()
loc.pca.plots <- list()

for (i in times){
  omic.time <-  miRNA[miRNA$Time_Point == i,]
  omic.pcr <- prcomp(omic.time[,51:(ncol(omic.time))], center = TRUE,scale. = TRUE)
  group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group") + ggtitle(i)
  loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
}


png("miRNA_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
dev.off()

png("miRNA_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
dev.off()


### adductomics ###


adductomics <- inner_join(Full_covar_add, adductomics, by = "XP_ID")


times <- levels(adductomics$Per)

group.pca.plots <- list()
loc.pca.plots <- list()

for (i in times){
  omic.time <-  adductomics[adductomics$Per == i,]
  omic.pcr <- prcomp(omic.time[,61:(ncol(omic.time))], center = TRUE,scale. = TRUE)
  group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group.x") + ggtitle(i)
  loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
}


png("adductomics_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
dev.off()

png("adductomics_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
dev.off()



################################################
################################################

### For significant OMICs

#### Metabolites#####

metabolites.sig.dat <- inner_join(Full_covar_add, metabolites.sig.dat, by = "XP_ID")

times <- levels(metabolites.sig.dat$Per)



group.pca.plots <- list()
loc.pca.plots <- list()

for (i in times){
  omic.time <-  metabolites.sig.dat[metabolites.sig.dat$Per == i,]
  omic.pcr <- prcomp(omic.time[,45:(ncol(omic.time))], center = TRUE,scale. = TRUE)
  group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group") + ggtitle(i)
  loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
}


png("Sig_metabolites_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
dev.off()

png("Sig_metabolites_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
dev.off()


### mRNA ###
# 
# mRNA.sig.dat <- inner_join(Full_covar_add, mRNA.sig.dat, by = "XP_ID")
# 
# times <- levels(mRNA.sig.dat$Time_Point)
# 
# group.pca.plots <- list()
# loc.pca.plots <- list()
# 
# 
# for (i in times){
#   omic.time <-  mRNA.sig.dat[mRNA.sig.dat$Time_Point == i,]
#   omic.pcr <- prcomp(omic.time[,59:(ncol(omic.time))], center = TRUE,scale. = TRUE)
#   group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group") + ggtitle(i)
#   loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
# }
# 
# 
# png("Sig_mRNA_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
# grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
# dev.off()
# 
# png("Sig_mRNA_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
# grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
# dev.off()


#### miRNA #####

miRNA.sig.dat <- inner_join(Full_covar_add, miRNA.sig.dat, by = "XP_ID")

times <- levels(miRNA.sig.dat$Per)

miRNA.sig.dat <- miRNA.sig.dat[-as.numeric(which(apply(miRNA.sig.dat, 2, var) == 0))] #removing columns with zero variance



group.pca.plots <- list()
loc.pca.plots <- list()

for (i in times){
  omic.time <-  miRNA.sig.dat[miRNA.sig.dat$Per == i,]
  omic.pcr <- prcomp(omic.time[,43:(ncol(omic.time))], center = TRUE,scale. = TRUE)
  group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group") + ggtitle(i)
  loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
}


png("Sig_miRNA_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
dev.off()

png("Sig_miRNA_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
dev.off()


### adductomics ###

# 
# adductomics.sig.dat <- inner_join(Full_covar_add, adductomics.sig.dat, by = "XP_ID")
# 
# 
# times <- levels(adductomics.sig.dat$Per)
# 
# group.pca.plots <- list()
# loc.pca.plots <- list()
# 
# for (i in times){
#   omic.time <-  adductomics.sig.dat[adductomics.sig.dat$Per == i,]
#   omic.pcr <- prcomp(omic.time[,61:(ncol(omic.time))], center = TRUE,scale. = TRUE)
#   group.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Group.x") + ggtitle(i)
#   loc.pca.plots[[i]] <- autoplot(omic.pcr, data = omic.time, colour = "Loc") + ggtitle(i)
# }
# 
# 
# png("Sig_adductomics_pca_loc_3timepoints.png",width = 6000, height = 2000, res = 300)
# grid.arrange(loc.pca.plots[[1]],loc.pca.plots[[2]],loc.pca.plots[[3]], ncol=3)
# dev.off()
# 
# png("Sig_adductomics_pca_group_3timepoints.png",width = 6000, height = 2000, res = 300)
# grid.arrange(group.pca.plots[[1]],group.pca.plots[[2]],group.pca.plots[[3]], ncol=3)
# dev.off()