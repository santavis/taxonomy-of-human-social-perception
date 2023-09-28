# Calculate significance of PCoA solutions for the generalization datasets independently and then compare the results

library(stringr)
library(corrplot)
library(lessR)

# Calculate the real PCoA solutions first for both datasets

# Load data
load("/path/clustering_analyses/analysis2.RData")
data_movie <- data
load("/path/clustering_analyses/analysis3.RData")
data_clips <- data

data_movie <- data_movie[,colnames(data_clips)]

# Distance measure
distance_clips <- 1-cor(data_clips) #Pearson
distance_movie <- 1-cor(data_movie) #Pearson

# PCoA with real data
fit_clips <- cmdscale(distance_clips,eig=TRUE,k=75)
fit_movie <- cmdscale(distance_movie,eig=TRUE,k=75)
loadings_clips <- fit_clips$points
loadings_movie <- fit_movie$points
weights_movie <- fit_movie$eig
weights_clips <- fit_clips$eig
var_exp_movie <- weights_movie[weights_movie>0]/sum(weights_movie[weights_movie>0])
var_exp_clips <- weights_clips[weights_clips>0]/sum(weights_clips[weights_clips>0])

##----------------------------------------------------------------------------------------------------------------------------------------------------
# How many significant components based on the clips 76 features?

path <- "/path/pcoa/permutations_generalization_primaryset"
batch <- 1000 # Batch size
rounds <- length(dir(path))/9 # How many batches have already been calculated
p <- batch*rounds

# Load null distributions
weigths_perm <- matrix(0,nrow=ncol(data_movie),ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("Movie: Loading null dist: ",i,"/",rounds,sep = ""))
  w <- read.csv(paste(path,"/pcoa_components_nulldist_",as.integer(i*batch),".csv",sep =""))
  weigths_perm[,((i-1)*batch+1):(i*batch)] <- as.matrix(w[,2:(batch+1)])
}

# Significance of components
s_component_clips <- c()
p_component_clips <- c()
for(j in seq(from=1, to=nrow(weigths_perm))){
  s_component_clips[j] <- sum((weigths_perm[j,] > weights_clips[j]),na.rm = TRUE)
  p_component_clips[j] <- s_component_clips[j]/p # This gives the exact p-value for each PC
}

##----------------------------------------------------------------------------------------------------------------------------------------------------
# How many significant components based on the movie 76 features?

path <- "/path/pcoa/permutations_generalization_validationset"
batch <- 1000 # Batch size
rounds <- length(dir(path))/9 # How many batches have already been calculated
p <- batch*rounds

# Load null distributions
weigths_perm <- matrix(0,nrow=ncol(data_movie),ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("Movie: Loading null dist: ",i,"/",rounds,sep = ""))
  w <- read.csv(paste(path,"/pcoa_components_nulldist_",as.integer(i*batch),".csv",sep =""))
  weigths_perm[,((i-1)*batch+1):(i*batch)] <- as.matrix(w[,2:(batch+1)])
  
}

# Significance of components
s_component_movie <- c()
p_component_movie <- c()
for(j in seq(from=1, to=nrow(weigths_perm))){
  s_component_movie[j] <- sum((weigths_perm[j,] > weights_movie[j]),na.rm = TRUE)
  p_component_movie[j] <- s_component_movie[j]/p # This gives the exact p-value for each PC
}

##---------------------------------------------------------------------------------------------------------------------------------------------------
# Compare and plot the correlations between the weights of the significant PCs.

# Take first 7 dimensions from clips and 10 dimensions from movie (significant components)
loadings7_clips <- loadings_clips[,1:7]
cats <- c("Clips PC1","Clips PC2","Clips PC3","Clips PC4","Clips PC5","Clips PC6","Clips PC7") 
colnames(loadings7_clips) <- cats

loadings10_movie <- loadings_movie[,1:10]
loadings10_movie[,1] <- loadings10_movie[,1]*-1 # Flip weight signs to avoid high negative correlations in the correlations matrix. The direction of the axis is insignificant
loadings10_movie[,5] <- loadings10_movie[,5]*-1 # Flip weight signs to avoid high negative correlations in the correlations matrix. The direction of the axis is insignificant
loadings10_movie[,9] <- loadings10_movie[,9]*-1 # Flip weight signs to avoid high negative correlations in the correlations matrix. The direction of the axis is insignificant
cats <- c("Movie PC1","Movie PC2","Movie PC3","Movie PC4","Movie PC5","Movie PC6","Movie PC7","Movie PC8","Movie PC9","Movie PC10")
colnames(loadings10_movie) <- cats

loadings <- cbind(loadings7_clips,loadings10_movie)
cormat <- corReorder(cor(loadings),order = "hclust")
loadings <- loadings[,colnames(cormat)]
cortest <- cor.mtest(loadings, conf.level = 0.95)
