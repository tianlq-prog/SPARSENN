## load positive data
load("ST001637 positive HMDB 5e-05 5 0.1 1 104 0.1 200 .rda")

pos_data <- g

## load negative data)
load("ST001637 negative HMDB 5e-05 5 0.1 1 104 0.1 200 .rda")

neg_data <- g


# delete some samples with no detail information
pos_data <- pos_data[,((!grepl('pool', colnames(pos_data), fixed = TRUE)&(!grepl('Mtdblk', colnames(pos_data), fixed = TRUE))))]
neg_data <- neg_data[,((!grepl('pool', colnames(neg_data), fixed = TRUE)&(!grepl('Mtdblk', colnames(neg_data), fixed = TRUE))))]


dim(pos_data) # 10085   484
dim(neg_data) # 6947  484

colnames(pos_data) <- gsub('pos', '', colnames(pos_data))
colnames(neg_data) <- gsub('neg', '', colnames(neg_data))

rownames(pos_data) <- paste('pos.',rownames(pos_data), sep = "")
rownames(neg_data) <- paste('neg.',rownames(neg_data), sep = "")

data <- rbind(pos_data, neg_data)
dim(data) #  17032   484


# about zero expression
rowsum <- rowSums(data==0)
hist(rowsum, breaks = 30)
summary(rowsum[rowsum>0])
sum(rowsum< 0.75 * dim(data)[2])/(dim(data)[1]) 


data <- data[which(rowsum<0.75 * dim(data)[2]),]

dim(data) # 16120   484

load("mouse_annotation.Rdata")

data_idx <- which(data[,1] %in% annotation[,9])
data_anno <- data[data_idx,] # data with annotation

dim(data_anno)  # 4209  484

annotation_1 <- annotation[which(annotation[,9] %in% data_anno[,1]),]
dim(annotation_1) # 13201     10

# the compound network
load('human.graph.degree.less.than.20.RData')
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
length(sel_compound) # 713

idx_g <- c()
for (i in sel_compound) {
  idx_g <- c(idx_g, which(all_compound==i))
}

adj_new <- as_adjacency_matrix(g)[idx_g,][,idx_g]

dim(adj_new) # 713 713

annotation_2 <- annotation_1[which(annotation_1$KEGGID %in% sel_compound),]
rownames(annotation_2) <- 1:nrow(annotation_2)  
dim(annotation_2) # 1113 10

data_anno_new <- data_anno[which(data_anno[,1] %in% annotation_2[,9]),]
dim(data_anno_new) # 671 484


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


dim(data_anno_new) # 671 484
dim(matching) # 671 713
dim(adj_new) # 713 713
dim(ion_matrix) # 671 713

write.csv(matching, "processed/feature_meta_matching.csv")
write.csv(data_anno_new, "processed/data.csv")
write.csv(as.matrix(adj_new), "processed/adj_mtx.csv")
write.csv(ion_matrix, "processed/ion_matching.csv")




