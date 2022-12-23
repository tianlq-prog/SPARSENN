rm(list = ls())
getwd()
load("annotation.RData")

load("raw_data/positive batchwise HMDB 12 0.5 2 0.333333333333333 0.5 0.5 new.Rdata")

pos_data <- g

## load negative data
load("raw_data/negative batchwise HMDB 12 0.5 2 0.333333333333333 0.5 0.5 new.Rdata")

neg_data <- g

#K <- rbind(K_pos, K_neg)

#dim(K) # 23206 10

dim(pos_data) # 4593 * 704
dim(neg_data) # 2733 * 704

colnames(pos_data) <- gsub('pos', '', colnames(pos_data))
colnames(neg_data) <- gsub('neg', '', colnames(neg_data))

rownames(pos_data) <- paste('pos.',rownames(pos_data), sep = "")
rownames(neg_data) <- paste('neg.',rownames(neg_data), sep = "")

data <- rbind(pos_data, neg_data)
dim(data) # 7326 * 704

# patient information
load("raw_data/polar_pos_batch_info.Rdata")
a$local_sample_id <- gsub('pos','', a$local_sample_id)
idx <- c()
for (n in colnames(data)) {idx <- c(idx, which(a$local_sample_id==n))}
batch <- a$batch[idx]
cov <- a$SARS.CoV.2.Positive[idx]
icu <- a$Admitted.to.the.ICU[idx]
day <- a$WU.day.of.presentation[idx]
Y_info <-data.frame(batch = batch, icu = icu, cov = cov, day = day) # information corresponding to Y
dim(Y_info) # 700 * 4

# about zero expression
rowsum <- rowSums(data==0)
hist(rowsum, breaks = 30)
summary(rowsum[rowsum>0])
sum(rowsum< 0.75 * 700)/(dim(data)[1]) 

data <- data[which(rowsum<0.75 * 700),]

dim(data) # 7326 * 700

# remove batch effect 
library(sva)
library(bladderbatch)
library(preprocessCore)
edata <- log(data[,5:dim(data)[2]]+1) # log transformation
combat_edata <- ComBat(dat = edata, batch = batch) #remove batch effect
combat_edata <- normalize.quantiles(combat_edata)
data[,5:dim(data)[2]] <- combat_edata

dim(data) # 7326 * 704

data_idx <- which(data[,1] %in% annotation[,9])
data_anno <- data[data_idx,] # data with annotation

dim(data_anno)  # 4308 * 704

annotation_1 <- annotation[which(annotation[,9] %in% data_anno[,1]),]
dim(annotation_1) # 23206 * 10


# the compound network
load("raw_data/human.graph.degree.less.than.20.RData")
library(igraph)
g <- as.undirected(g) # change to undirected
g_name <- V(g)$name
r_name <- g_name[substring(g_name,1,1)=='R']

add_edge_func <- function(g, name){
  g_name <- V(g)$name
  idx <- which(g_name==name)
  nb <- unique(neighbors(g, idx)$name)
  if (length(nb)>1){
    for (i in 1:(length(nb)-1)){
      for (j in (i+1):length(nb)){
        g <- add_edges(g, c(which(V(g)$name==nb[i]), which(V(g)$name==nb[j])), color = 'red')
      }
    }
  }
  g <- delete.vertices(g, which(V(g)$name==name))
  print(length(V(g)$name))
  return(g)
}

for (r in r_name){g <- add_edge_func(g, r)}


all_compound <- V(g)$name
sel_compound <- intersect(all_compound, unique(annotation_1$KEGGID))
length(sel_compound) # 913

idx_g <- c()
for (i in sel_compound) {
  idx_g <- c(idx_g, which(all_compound==i))
}

adj_new <- as_adjacency_matrix(g)[idx_g,][,idx_g]

dim(adj_new) # 913 * 913

annotation_2 <- annotation_1[which(annotation_1$KEGGID %in% sel_compound),]
rownames(annotation_2) <- 1:nrow(annotation_2)  
dim(annotation_2) # 2851 * 10

data_anno_new <- data_anno[which(data_anno[,1] %in% annotation_2[,9]),]
dim(data_anno_new) # 1351 * 704


# obtain the feature-meta-matrix
feature_ls <- data_anno_new[,1]
matching <- matrix(0, nrow = length(feature_ls), ncol = length(sel_compound))
ion_matrix <- matrix('0', nrow = length(feature_ls), ncol = length(sel_compound))
count = 1
for (i in sel_compound){
  sel_feature <- annotation_2[which(annotation_2$KEGGID== i),9]
  idx <- c()
  for (j in sel_feature) {
    idx <- c(idx, which(feature_ls==j))
  }
  matching[idx, count] <- 1
  ion_matrix[idx, count] <- as.character(annotation_2$Adduct[which(annotation_2$KEGGID==i)])
  count <- count + 1
}

colnames(matching) <- sel_compound
rownames(matching) <- rownames(data_anno_new)

colnames(ion_matrix) <- sel_compound
rownames(ion_matrix) <- rownames(data_anno_new)


dim(data_anno_new) # 1351 * 704
dim(matching) # 1351 * 913
dim(adj_new) # 913 * 913
dim(ion_matrix) # 1351 * 913



write.csv(matching, "processed/feature_meta_matching.csv")
write.csv(ion_matrix, "processed/ion_matching.csv")
write.csv(data_anno_new, "processed/data.csv")
write.csv(as.matrix(adj_new), "processed/adj_mtx.csv")
write.csv(Y_info, "processed/Y.csv")

