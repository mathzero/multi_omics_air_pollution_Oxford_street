work.dir <- "~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/"
setwd(work.dir)

exposure_significant<-read.csv('Results/phase_1/exposure_association_results.csv')


source("http://www.bioconductor.org/biocLite.R")
class(biocLite)
library(limma)
library(VennDiagram)
library(grid)
library(futile.logger)
library(gridBase)
library(lattice)
library(gridExtra)

attach(exposure_significant)


############# Loop for four OMICs ###################

omics <- c("mRNA","miRNA","metabolite","adductomic")

venn.list <- list()

for (i in omics){
  PM10<-subset(exposure_significant, OMIC_level==i & TRAP=="PM10")[,'OMIC.ID']
  NO2<-subset(exposure_significant, OMIC_level==i & TRAP=="NO2")[,'OMIC.ID']
  PM25<-subset(exposure_significant, OMIC_level==i & TRAP=="PM25")[,'OMIC.ID']
  CBLK<-subset(exposure_significant, OMIC_level==i & TRAP=="CBLK")[,'OMIC.ID']
  PCNT<-subset(exposure_significant, OMIC_level==i & TRAP=="PCNT")[,'OMIC.ID']
  
  mysets= list(PM10=PM10, NO2=NO2, PM25=PM25,CBLK=CBLK, PCNT=PCNT) 
  
  colors = c("tomato", "forestgreen", "skyblue", "gold",
             "navy")
  
  names = c("PM10", "NO2", "PM25", "CBLK", "PCNT")
  
  venn.list[[i]] <- venn.diagram(mysets, filename = NULL,
                 fill = colors, category.names = names, 
                 main = i,
                 lwd=0,
                 alpha =0.4,
                 fontfamily = "sans",
                 cat.fontfamily = "sans",
                 main.fontfamily = "sans",
                 main.cex = 1.6,
                 cat.just = list(c(0.5, 1), c(-0.5, -5), c(0, 0),
                                 c(0.5, 0), c(1, -4)))

}

###### Creating combined plot #######

png("venn_diagrams_combined.png", width = 3000, height = 3000, res = 300)
grid.arrange(grobTree(venn.list[[1]]), grobTree(venn.list[[2]]),
             grobTree(venn.list[[3]]), grobTree(venn.list[[4]]), 
             ncol=2, padding = unit(5, "line"))
dev.off()

