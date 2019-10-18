work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)

library(mixOmics)


metabolites_significant_denoised <- readRDS("Data/phase_2/denoised/metabolites_significant_denoised.rds")
adductomics_significant_denoised <- readRDS("Data/phase_2/denoised/adductomics_significant_denoised.rds")
miRNA_significant_denoised <- readRDS("Data/phase_2/denoised/miRNA_significant_denoised.rds")
mRNA_significant_denoised <- readRDS("Data/phase_2/denoised/mRNA_significant_denoised.rds")
mRNA_denoised <- readRDS("Data/phase_2/denoised/mRNA_denoised.rds")
adductomics_denoised <- readRDS("Data/phase_2/denoised/adductomics_denoised.rds")
metabolites_denoised <- readRDS("Data/phase_2/denoised/metabolites_denoised.rds")
miRNA_denoised <- readRDS("Data/phase_2/denoised/miRNA_denoised.rds")

dat <- cbind(mRNA_significant_denoised, 
             miRNA_significant_denoised[,2:ncol(miRNA_significant_denoised)], 
             metabolites_significant_denoised[,2:ncol(metabolites_significant_denoised)],
             adductomics_significant_denoised[,2:ncol(adductomics_significant_denoised)])

dat <- dat[dat$id != 31,]

table(dat$id)
#For future ref:

# dat 1:4  =  loc, XP_ID, id, per
# dat 5    =  mRNA
# dat 6:11 =  miRNA
# dat 12:160   =  metabolites
# dat 160:162  =  adductomics



### Tuning ###
tune = tune.spls(as.matrix(dat[,6:11]), as.matrix(dat[,12:160]), ncomp=2, test.keepX = c(1:20), measure = "MSE", 
                 nrepeat=3,progressBar = TRUE)

tune$choice.ncomp
tune$choice.keepX

# plot the results
plot(tune)
tune$choice.ncomp



### Building model ###
spls.mod <- spls(as.matrix(dat[,6:11]),
                 as.matrix(dat[,5]),
     ncomp = 1,
     mode ="classic",
     scale = TRUE,
     tol = 1e-06,
     max.iter = 100,
     near.zero.var = FALSE,
     logratio="none", 
     multilevel=as.matrix(dat[,3]),
     all.outputs = TRUE)

plotLoadings(spls.mod,contrib = 'max', comp = 1, method = 'median', title = "sPLS mRNA on miRNA")


### Building model ###

spls.mod <- spls(as.matrix(dat[,12:160]),
                 as.matrix(dat[,6:11]),
                 ncomp = 1,
                 mode ="classic",
                 scale = TRUE,
                 tol = 1e-06,
                 max.iter = 100,
                 near.zero.var = FALSE,
                 logratio="none", 
                 multilevel=as.matrix(dat[,3]),
                 all.outputs = TRUE,
                 keepX = c(10,1), keepY = c(5,5))

plotLoadings(spls.mod,contrib = 'max', comp = 1, method = 'median', title = "sPLS mRNA on miRNA")

val <- perf(spls.mod, criterion = c("Q2","R2"), validation = "loo")
plot(val$Q2)




