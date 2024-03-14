# Permutation test for the significance of the PCoA components of the primary movie clip dataset (Study 2).
# The null distributions are calculated previously.

# Concordance analysis of Study 2. PCoA component loadings are calculated for each hierarchical cluster to estimate the relationship between the two methods.

# Severi Santavirta 26.2.2024

library(stringr)

# P-value function (one-sided test)
pvalue <- function(value,nulldist){
  if(value>0){
    pval <- sum(nulldist > value)/length(nulldist)
  } else {
    pval <- sum(nulldist < value)/length(nulldist)
  }
  return(pval)
}

# Permutation folder
path <- "/path/permutations"
batch <- 1000 # Batch size
rounds <- length(dir(path))/9 # How many batches have already been calculated
p <- batch*rounds

# Load data
tbl1 <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_1.csv")
tbl2 <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_2.csv")
tbl3 <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_3.csv")
tbl4 <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_4.csv")
tbl5 <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_5.csv")
tbl6 <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_6.csv")
data <- rbind(tbl1,tbl2,tbl3,tbl4,tbl5,tbl6)

dims <- colnames(data)
dim_names <- str_replace_all(dims,"_"," ")

# Select Pearson correlation distance as distance measure
distance <- 1-cor(data)

# Run PCoA on real data
fit <- cmdscale(distance,eig=TRUE,k=135)
loadings <- fit$points
weights <- fit$eig
var_exp <- weights[weights>0]/sum(weights[weights>0])

# Load null distributions
weigths_perm <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm1 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm2 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm3 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm4 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm5 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm6 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm7 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm8 <- matrix(0,nrow=136,ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("Loading null dist: ",i,"/",rounds,sep = ""))
  w <- read.csv(paste(path,"/pcoa_components_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc1 <- read.csv(paste(path,"/pcoa_loadings_pc1_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc2 <- read.csv(paste(path,"/pcoa_loadings_pc2_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc3 <- read.csv(paste(path,"/pcoa_loadings_pc3_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc4 <- read.csv(paste(path,"/pcoa_loadings_pc4_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc5 <- read.csv(paste(path,"/pcoa_loadings_pc5_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc6 <- read.csv(paste(path,"/pcoa_loadings_pc6_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc7 <- read.csv(paste(path,"/pcoa_loadings_pc7_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc8 <- read.csv(paste(path,"/pcoa_loadings_pc8_nulldist_",as.integer(i*batch),".csv",sep =""))
  
  weigths_perm[,((i-1)*batch+1):(i*batch)] <- as.matrix(w[,2:(batch+1)])
  loadings_perm1[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc1[,2:(batch+1)])
  loadings_perm2[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc2[,2:(batch+1)])
  loadings_perm3[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc3[,2:(batch+1)])
  loadings_perm4[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc4[,2:(batch+1)])
  loadings_perm5[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc5[,2:(batch+1)])
  loadings_perm6[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc6[,2:(batch+1)])
  loadings_perm7[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc7[,2:(batch+1)])
  loadings_perm8[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc8[,2:(batch+1)])
}

# Significance of components
s_component <- c()
p_component <- c()
for(j in seq(from=1, to=135)){
  s_component[j] <- sum((weigths_perm[j,] > weights[j]),na.rm = TRUE)
  p_component[j] <- s_component[j]/p
}

# Significance of loadings
p_loadings <- matrix(0,nrow = 136, ncol = 8)
for(j in seq(from=1, to=136)){
  
  if(loadings[j,1] > 0) {
    p_loadings[j,1] <- sum((loadings_perm1[j,] > loadings[j,1]), na.rm = TRUE)/p
  } else {
    p_loadings[j,1] <- sum((loadings_perm1[j,] < loadings[j,1]), na.rm = TRUE)/p
  }
  
  if(loadings[j,2] > 0) {
    p_loadings[j,2] <- sum((loadings_perm2[j,] > loadings[j,2]), na.rm = TRUE)/p
  } else {
    p_loadings[j,2] <- sum((loadings_perm2[j,] < loadings[j,2]), na.rm = TRUE)/p
  }
  
  if(loadings[j,3] > 0) {
    p_loadings[j,3] <- sum((loadings_perm3[j,] > loadings[j,3]), na.rm = TRUE)/p
  } else {
    p_loadings[j,3] <- sum((loadings_perm3[j,] < loadings[j,3]), na.rm = TRUE)/p
  }
  
  if(loadings[j,4] > 0) {
    p_loadings[j,4] <- sum((loadings_perm4[j,] > loadings[j,4]), na.rm = TRUE)/p
  } else {
    p_loadings[j,4] <- sum((loadings_perm4[j,] < loadings[j,4]), na.rm = TRUE)/p
  }
  
  if(loadings[j,5] > 0) {
    p_loadings[j,5] <- sum((loadings_perm5[j,] > loadings[j,5]), na.rm = TRUE)/p
  } else {
    p_loadings[j,5] <- sum((loadings_perm5[j,] < loadings[j,5]), na.rm = TRUE)/p
  }
  
  if(loadings[j,6] > 0) {
    p_loadings[j,6] <- sum((loadings_perm6[j,] > loadings[j,6]), na.rm = TRUE)/p
  } else {
    p_loadings[j,6] <- sum((loadings_perm6[j,] < loadings[j,6]), na.rm = TRUE)/p
  }
  
  if(loadings[j,7] > 0) {
    p_loadings[j,7] <- sum((loadings_perm7[j,] > loadings[j,7]), na.rm = TRUE)/p
  } else {
    p_loadings[j,7] <- sum((loadings_perm7[j,] < loadings[j,7]), na.rm = TRUE)/p
  }
  
  if(loadings[j,8] > 0) {
    p_loadings[j,8] <- sum((loadings_perm8[j,] > loadings[j,8]), na.rm = TRUE)/p
  } else {
    p_loadings[j,8] <- sum((loadings_perm8[j,] < loadings[j,8]), na.rm = TRUE)/p
  }
  
}

rownames(p_loadings) <- dim_names
sign001_loadings <- (p_loadings < 0.01)
signbonf_loadings <- (p_loadings < 0.05/136) # Correct for number of tested dimensions

#write.csv(p_component,"/path/pcoa_components_pval.csv")
#write.csv(p_loadings,"/path/pcoa_loadings_pval.csv")
#write.csv(sign001_loadings,"/path/pcoa_loadings_pval_sign001.csv")
#write.csv(signbonf_loadings,"/path/pcoa_loadings_pval_signbonf.csv")

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Significance of average loadings within each hierarchical cluster (Study 2: Concordance analysis)

library(lessR)

load("/path/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

distance <- 1-cor(data) #Pearson

fit <- cmdscale(distance,eig=TRUE,k=135)
loadings <- fit$points

#Cluster membership by cluster
cc <- c(rep(1,136))
cc[2] <- 2
cc[3:5] <- 3
cc[6:7] <- c(4,5)
cc[8:11] <- 6
cc[12] <- 7
cc[13:16] <- 8
cc[17:26] <- 9
cc[27:32] <- 10
cc[33:50] <- 11
cc[51] <- 12
cc[52:57] <- 13
cc[58:59] <- 14
cc[60] <- 15
cc[61:66] <- 16
cc[67] <- 17
cc[68:70] <- 18
cc[71:72] <- 19
cc[73:79] <- 20
cc[80:83] <- 21
cc[84:101] <- 22
cc[102:115] <- 23
cc[116] <- 24
cc[117:121] <- 25
cc[122:123] <- 26
cc[124:128] <- 27
cc[129:130] <- 28
cc[131] <- 29
cc[132:136] <- 30

names <- colnames(data)
df <- data.frame(names,loadings)
df$group <- factor(cc,levels = unique(cc))

# No clusters
idx <- (df$group==1 | df$group==2 | df$group==4 | df$group==5 | df$group==7 | df$group==12 | df$group==15 | df$group==17 | df$group==24 | df$group==29)
df1 <- df[!idx,]

# Load null distributions
path <- "/path/permutations"
batch <- 1000 # Batch size
rounds <- length(dir(path))/9 # How many batches have already been calculated
p <- batch*rounds

# Load null distributions

loadings_perm1 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm2 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm3 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm4 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm5 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm6 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm7 <- matrix(0,nrow=136,ncol=batch*rounds)
loadings_perm8 <- matrix(0,nrow=136,ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("Loading null dist: ",i,"/",rounds,sep = ""))
  pc1 <- read.csv(paste(path,"/pcoa_loadings_pc1_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc2 <- read.csv(paste(path,"/pcoa_loadings_pc2_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc3 <- read.csv(paste(path,"/pcoa_loadings_pc3_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc4 <- read.csv(paste(path,"/pcoa_loadings_pc4_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc5 <- read.csv(paste(path,"/pcoa_loadings_pc5_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc6 <- read.csv(paste(path,"/pcoa_loadings_pc6_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc7 <- read.csv(paste(path,"/pcoa_loadings_pc7_nulldist_",as.integer(i*batch),".csv",sep =""))
  pc8 <- read.csv(paste(path,"/pcoa_loadings_pc8_nulldist_",as.integer(i*batch),".csv",sep =""))
  
  loadings_perm1[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc1[,2:(batch+1)])
  loadings_perm2[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc2[,2:(batch+1)])
  loadings_perm3[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc3[,2:(batch+1)])
  loadings_perm4[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc4[,2:(batch+1)])
  loadings_perm5[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc5[,2:(batch+1)])
  loadings_perm6[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc6[,2:(batch+1)])
  loadings_perm7[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc7[,2:(batch+1)])
  loadings_perm8[,((i-1)*batch+1):(i*batch)] <- as.matrix(pc8[,2:(batch+1)])
}

# Calculate mean loading and its significance for each dimension and each cluster
#   Significance is calculated by creating the null distributions for average ratings based on the previously calculated null_distributions for each feature
dim_order <- pc1$X
dim_order[92] <- "Exerting self-control" 
idx <- unique(df1$group)
p_cluster <- matrix(1,ncol=length(idx),nrow = 8)
for(i in seq(from=1,to=length(idx))){
  data_cluster <- df1[df1$group==idx[i],]
  avg_cluster <- colMeans(data_cluster[,2:9]) #True average
  
  idx_col <- which(dim_order %in% data_cluster$names)
  if(length(idx_col)!=nrow(data_cluster)) {
    print(paste("Problem!",i))
  }
  
  # Null distributions of the average loading
  avg1_null <- colMeans(loadings_perm1[idx_col,])
  avg2_null <- colMeans(loadings_perm2[idx_col,])
  avg3_null <- colMeans(loadings_perm3[idx_col,])
  avg4_null <- colMeans(loadings_perm4[idx_col,])
  avg5_null <- colMeans(loadings_perm5[idx_col,])
  avg6_null <- colMeans(loadings_perm6[idx_col,])
  avg7_null <- colMeans(loadings_perm7[idx_col,])
  avg8_null <- colMeans(loadings_perm8[idx_col,])
  
  #P-value
  p_cluster[1,i] <- pvalue(avg_cluster[1],avg1_null)
  p_cluster[2,i] <- pvalue(avg_cluster[2],avg2_null)
  p_cluster[3,i] <- pvalue(avg_cluster[3],avg3_null)
  p_cluster[4,i] <- pvalue(avg_cluster[4],avg4_null)
  p_cluster[5,i] <- pvalue(avg_cluster[5],avg5_null)
  p_cluster[6,i] <- pvalue(avg_cluster[6],avg6_null)
  p_cluster[7,i] <- pvalue(avg_cluster[7],avg7_null)
  p_cluster[8,i] <- pvalue(avg_cluster[8],avg8_null)
  
}

sign001_clusters <- (p_cluster < 0.01)
signbonf_clusters <- (p_cluster < 0.05/20) # Correct for number of tested clusters

p_cluster <- as.data.frame(p_cluster)
sign001_clusters <- as.data.frame(sign001_clusters)
signbonf_clusters <- as.data.frame(signbonf_clusters)

clusters <- c("Introversion","Feeding","Self-control","Achievement","Emotional affection","Pleasant feelings & prosociality","Extraversion & playfulness","Masculinity","Social engagement","Physical discomfort","Gesturing","Emotional expression","Motivation","Antisocial behaviour","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical interaction")
colnames(p_cluster) <- clusters
colnames(sign001_clusters) <- clusters
colnames(signbonf_clusters) <- clusters

#Save for later use 
#write.csv(p_cluster,"/path/pcoa_cluster_loadings_pval.csv")
#write.csv(sign001_clusters,"/path/pcoa_cluster_loadings_pval_sing001.csv")
#write.csv(signbonf_clusters,"/path/pcoa_cluster_loadings_pval_signbonf.csv")
