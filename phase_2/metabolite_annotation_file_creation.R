metabolites <- read.csv("~/Google Drive/Imperial/2 Translational data science/TDS project/Project6/Results/phase_1/results_metabolomics/metabolites_all_pvals.csv")


metabs <- metabolites[,2:3]

colnames(metabs) <- c("m.z", "p.value")


metabs$t.score <- (1/metabs$p.value)^0.2



metabs$m.z <- gsub("X","",metabs$m.z)
metabs$m.z <- sub("^(.*)[.].*", "\\1", metabs$m.z)
metabs$m.z <- sub("^(.*)[.].*", "\\1", metabs$m.z)
metabs$m.z <- substr(metabs$m.z, 1, nchar(metabs$m.z)-2) 

metabs$p.value <- round(metabs$p.value,5)
metabs$t.score <- round(metabs$t.score,5)

write.table(metabs[-c(1056,1057,3308,3309,3351,3352),], "metab_annotation.txt",row.names=FALSE)



metabs[metabs$m.z == "558.7933",]
