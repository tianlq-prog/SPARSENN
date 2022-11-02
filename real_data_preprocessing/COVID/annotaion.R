rm(list = ls())
getwd()

library(xMSannotator)
load("positive batchwise HMDB 12 0.5 2 0.333333333333333 0.5 0.5 new.Rdata")
g_pos <- g

colnames(g_pos)[1] <- "m/z"
colnames(g_pos)[2] <- "time"

# positive mode
K_positive<-KEGG.Annotation(g_pos[,1:4], max.mz.diff = 5, num_nodes = 2,
                            queryadductlist = c("M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H"),  
                            mode = "pos",
                            outloc = getwd()) #, mode = "pos", outloc=getwd()

K_positive$mz <- as.numeric(as.character(K_positive$mz))
K_positive$time <- as.numeric(as.character(K_positive$time))
K_positive$MonoisotopicMass <- as.numeric(as.character(K_positive$MonoisotopicMass))
K_positive$AdductMass <- as.numeric(as.character(K_positive$AdductMass))

K_positive <- K_positive[!duplicated(K_positive),] # delete deuplicated part

# negative mode
load("negative batchwise HMDB 12 0.5 2 0.333333333333333 0.5 0.5 new.Rdata")
g_neg <- g

colnames(g_neg)[1] <- "m/z"
colnames(g_neg)[2] <- "time"

K_negative<-KEGG.Annotation(g_neg[,1:4], max.mz.diff = 5, num_nodes = 2,queryadductlist =
                              c("M-H", "M-2H", "M-2H+Na", "M-2H+K", "M-2H+NH4", 
                                "M-H2O-H", "M-H+Cl", "M+Cl", "M+2Cl"), 
                            mode = "neg", outloc=getwd())

K_negative$mz <- as.numeric(as.character(K_negative$mz))
K_negative$time <- as.numeric(as.character(K_negative$time))
K_negative$MonoisotopicMass <- as.numeric(as.character(K_negative$MonoisotopicMass))
K_negative$AdductMass <- as.numeric(as.character(K_negative$AdductMass))

K_negative <- K_negative[!duplicated(K_negative),] # delete deuplicated part

# combine two mode annotation
annotation <- rbind(K_positive, K_negative)
dim(annotation) # 23206    10

save(annotation, file = "annotation.Rdata")









