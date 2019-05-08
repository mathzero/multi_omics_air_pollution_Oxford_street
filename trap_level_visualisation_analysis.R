rm(list=ls())
work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)

library(dplyr)
library(ggplot2)
library(gridExtra)

################################################
################################################

omics <- colnames(dat.sig[3:162])


dat.sig <- read.csv("multiomic_significant_vars.csv")
dat.after1 <- dat.sig[dat.sig$Per == "after_1",]

plot.list <- list()



for (i in 1:160){
  omicvar = colnames(dat.after1[(i+2)]) 
  plot.list[[i]] <- 
    ggplot(dat.after1, aes_string(x = "Loc", y = omicvar)) +   
    geom_boxplot(aes(fill = Loc), alpha = 0.5) +
    geom_line(aes(group = id),
              alpha = 0.2, colour = "red") +
      geom_jitter(alpha=0.5, aes(color=Loc),
                  position = position_jitter(width = 0.1))+
      ggtitle(paste0("Expression: ",omicvar))
} 



################################################
################################################

#plotting and saving

#miRNAs
png("miRNA_expression_bySite.png",width = 3000, height = 5000, res = 300)
grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],ncol = 2)
dev.off()

#mRNAs
png("mRNA_expression_bySite.png",width = 3000, height = 3000, res = 300)
grid.arrange(plot.list[[7]],plot.list[[8]],ncol = 2)
dev.off()

#metabolites1
metab.list <- plot.list[c(9:85)]
png("metabolite_expression_bySite1.png",width = 8000, height = 14000, res = 300)
do.call("grid.arrange", c(metab.list, ncol=6))
dev.off()

#metabolites2
metab.list <- plot.list[c(86:157)]
png("metabolite_expression_bySite2.png",width = 8000, height = 14000, res = 300)
do.call("grid.arrange", c(metab.list, ncol=6))
dev.off()


#adductomics
png("adductomics_expression_bySite.png",width = 8000, height = 3000, res = 300)
grid.arrange(plot.list[[158]],plot.list[[159]],ncol = 2)
dev.off()



################################################
################################################

Full_covar_add$Per <- factor(Full_covar_add$Per, levels = c("before","after_1","after_2"))

g1 <- ggplot(Full_covar_add, aes(NO2_annual,NO2_r, col = Loc)) + 
  geom_point(size = 1, alpha = 0.6) + 
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, max(Full_covar_add$NO2_annual)) +
  facet_wrap(~Per) + 
  ggtitle("NO2 levels relative to baseline, at different time points")

g2 <- ggplot(Full_covar_add, aes(PM10_annual,PM10_r, col = Loc)) + 
  geom_point(size = 1, alpha = 0.6) + 
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, max(Full_covar_add$PM10_annual)) +
  facet_wrap(~Per) + 
  ggtitle("PM10 levels relative to baseline, at different time points")

g3 <- ggplot(Full_covar_add, aes(PM25_annual,PM25_r, col = Loc)) + 
  geom_point(size = 1, alpha = 0.6) + 
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, max(Full_covar_add$PM25_annual)) +
  facet_wrap(~Per) + 
  ggtitle("PM25 levels relative to baseline, at different time points")

g4 <- ggplot(Full_covar_add, aes(PCNT_annual,PCNT_r, col = Loc)) + 
  geom_point(size = 1, alpha = 0.6) + 
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, max(Full_covar_add$PCNT_annual)) +
  facet_wrap(~Per) + 
  ggtitle("UFP levels relative to baseline, at different time points")

g5 <- ggplot(Full_covar_add, aes(CBLK_annual,CBLK_r, col = Loc)) + 
  geom_point(size = 1, alpha = 0.6) + 
  geom_abline(intercept = 0, slope = 1, alpha = 0.3) +
  xlim(0, max(Full_covar_add$CBLK_annual)) +
  facet_wrap(~Per) + 
  ggtitle("Black Carbon levels relative to baseline, at different time points")


png("exposure_by_time.png",width = 3000, height = 6000, res = 300)
grid.arrange(g1,g2,g3,g4,g5, ncol=1)
dev.off()
