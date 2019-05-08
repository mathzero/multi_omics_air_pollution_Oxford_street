rm(list=ls())
work_dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project"

#"~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
#"/rds/general/user/mw418/home/TDS_Project"

setwd(work_dir)

library(dplyr)

################################################
################################################



### miRNA Manhattan plot ####


#Reading in miRNA chromosome info
mRNAchrom <- read.csv("~/Google Drive/Imperial/2 Translational data science/TDS project/Misc/group-476.csv", header=FALSE)
mRNAchrom <- mRNAchrom[3:nrow(mRNAchrom),]
colnames(mRNAchrom) <- c("HGNC_ID_gene","approved_symbol","approved_name","prev_symbols","synonym","chromosome")

#reading in miRNA full pvals #
miRNA_all_pvals <- read.csv("/Results/phase_1/results_miRNA/05_sensitivity/miRNA_all_pvals.csv",header=FALSE)

#manipulating format to join dfs #
miRNA_all_pvals$X <- substr(miRNA_all_pvals$OMICid,9,20)
miRNA_all_pvals$X <-tools::file_path_sans_ext(miRNA_all_pvals$X)
miRNA_all_pvals$X <- gsub("\\.", "-", miRNA_all_pvals$X)
miRNA_all_pvals$approved_symbol <- paste("MIR",miRNA_all_pvals$X,sep="")
miRNA_all_pvals$approved_symbol[1:9] <- gsub("MIR","MIRLET",miRNA_all_pvals$approved_symbol[1:9])
miRNA_all_pvals$approved_symbol <- toupper(miRNA_all_pvals$approved_symbol)

################################################
################################################

#joining dfs #

miRNA.manhat <- left_join(miRNA_all_pvals, mRNAchrom, by = "approved_symbol")
miRNA.manhat$chr <-   gsub( "[.].*$", "", miRNA.manhat$chromosome)
miRNA.manhat$chr <- gsub('.{3}$', '', miRNA.manhat$chr)
miRNA.manhat$chr[miRNA.manhat$chr == "X"]  <- 23
miRNA.manhat$chr <- as.numeric(miRNA.manhat$chr)
miRNA.manhat$bp <- as.numeric(rownames(miRNA.manhat))


# Load the library
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
manhattan(miRNA.manhat[!is.na(miRNA.manhat$chr),], chr="chr", bp="bp", snp="OMICid", p="NO2_pval",annotatePval = 0.05/nrow(miRNA.manhat))

cmplot1 <- miRNA.manhat[!is.na(miRNA.manhat$chr),c(2,14,15,3:7)]


library("CMplot")      #install.packages("CMplot")
CMplot(cmplot1, plot.type="m", multracks=TRUE, threshold=c(0.05,0.05/nrow(cmplot1)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=NULL, signal.col=c("red","green"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


