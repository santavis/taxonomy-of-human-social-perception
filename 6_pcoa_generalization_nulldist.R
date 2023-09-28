# Generalization analysis for the principal component solutions between the primary and validation datasets

library(stringr)

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

dim_names <- colnames(data_clips)

# Run permutations on randomized data to estimate chance distributions (first for movie)

p <- 1000000 # How many permutations
batch <- 1000
n <- 0 # set own seed for each sample
k <- 0 # Count permutations
weights_perm <- matrix(0,nrow = ncol(data_movie), ncol = batch)

for(i in seq(from=1,to=p)){
  k <- k+1
  if(i%%100 == 0){
    print(paste("Movie: permutation: ",i,"/",p))
  }
  data_perm <- data_movie
  for(j in seq(from=1,to=ncol(data_perm))){ # Create boostrapped sample
    n <- n+1
    set.seed(n)
    data_perm[,j] <- sample(data_perm[,j],replace = FALSE) # Selects all the data in a random order (based on the seed number, which is different for each column and each permutation)
  }
  distance_perm <- 1-cor(data_perm)
  fit_perm <- cmdscale(distance_perm,eig=TRUE,k=ncol(data_perm)-1)
  
  # Store weights to test the components significance
  weights_perm[,k] <- fit_perm$eig
  
  if(i%%batch == 0){ # Save null distributions for every batch
    
    rownames(weights_perm) <- dim_names
    write.csv(weights_perm,paste("/path/pcoa/permutations_generalization_validationset/pcoa_components_nulldist_",i,".csv",sep = ""))
    k<-0
  }
}

# Run permutations on randomized data to estimate chance distributions (then for clips)

p <- 1000000 # How many permutations
batch <- 1000
n <- 0 # set own seed for each sample
k <- 0 # Count permutations
weights_perm <- matrix(0,nrow = ncol(data_clips), ncol = batch)

for(i in seq(from=1,to=p)){
  k <- k+1
  if(i%%100 == 0){
    print(paste("clips: permutation: ",i,"/",p))
  }
  data_perm <- data_clips
  for(j in seq(from=1,to=ncol(data_perm))){ # Create boostrapped sample
    n <- n+1
    set.seed(n)
    data_perm[,j] <- sample(data_perm[,j],replace = FALSE) # Selects all the data in a random order (based on the seed number, which is different for each column and each permutation)
  }
  distance_perm <- 1-cor(data_perm)
  fit_perm <- cmdscale(distance_perm,eig=TRUE,k=ncol(data_perm)-1)
  
  # Store weights to test the components significance
  weights_perm[,k] <- fit_perm$eig
  
  if(i%%batch == 0){ # Save null distributions for every batch
    
    rownames(weights_perm) <- dim_names
    write.csv(weights_perm,paste("/path/pcoa/permutations_generalization_primaryset/pcoa_components_nulldist_",i,".csv",sep = ""))
    k<-0
  }
}

