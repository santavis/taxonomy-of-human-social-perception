# Principal coordinate analysis (PCoA) of social perceptual features with a permutation test for identifying significant components for the full movie generalization study (Study 3)

# Severi Santavirta 26.2.2024


##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate the null distributions for component weights for both the movie clip data (76 features) and for the full movie data (76 features)

library(stringr)

# Load data
load("/path/analysis2.RData")
data_movie <- data
load("/path/analysis3.RData")
data_clip <- data

data_movie <- data_movie[,colnames(data_clip)]

# Distance measure
distance_clip <- 1-cor(data_clip) #Pearson
distance_movie <- 1-cor(data_movie) #Pearson

# PCoA with real data
fit_clip <- cmdscale(distance_clip,eig=TRUE,k=75)
fit_movie <- cmdscale(distance_movie,eig=TRUE,k=75)
loadings_clip <- fit_clip$points
loadings_movie <- fit_movie$points

dim_names <- colnames(data_clip)

# Run permutations on randomized data to estimate chance distributions for full movie

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
    #write.csv(weights_perm,paste("/path/permutations_generalization_movie/pcoa_components_nulldist_",i,".csv",sep = ""))
    k<-0
  }
}

# Run permutations on randomized data to estimate chance distributions for movie clips

p <- 1000000 # How many permutations
batch <- 1000
n <- 0 # set own seed for each sample
k <- 0 # Count permutations
weights_perm <- matrix(0,nrow = ncol(data_clip), ncol = batch)

for(i in seq(from=1,to=p)){
  k <- k+1
  if(i%%100 == 0){
    print(paste("Clips: permutation: ",i,"/",p))
  }
  data_perm <- data_clip
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
    #write.csv(weights_perm,paste("/path/permutations_generalization_clips/pcoa_components_nulldist_",i,".csv",sep = ""))
    k<-0
  }
}


rm(list = ls())


##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate how many significant components movie clip data (76 features) and full movie data (76 features) have. (Study 3)
# This information is used when comparing the component loadings across datasets (see figures script)

# Severi Santavirta 26.2.2024

library(stringr)
library(corrplot)
library(lessR)

# Calculate the real PCoA solutions first for both datasets

# Load data
load("/path/analysis2.RData")
data_movie <- data
load("/path/analysis3.RData")
data_clip <- data

data_movie <- data_movie[,colnames(data_clip)]

# Distance measure
distance_clip <- 1-cor(data_clip) #Pearson
distance_movie <- 1-cor(data_movie) #Pearson

# PCoA with real data
fit_clip <- cmdscale(distance_clip,eig=TRUE,k=75)
fit_movie <- cmdscale(distance_movie,eig=TRUE,k=75)
loadings_clip <- fit_clip$points
loadings_movie <- fit_movie$points
weights_movie <- fit_movie$eig
weights_clip <- fit_clip$eig

##----------------------------------------------------------------------------------------------------------------------------------------------------
# How many significant components based on the 76 common features in movie clip dataset?

path <- "/path/permutations_generalization_clips"
batch <- 1000 # Batch size
rounds <- length(dir(path))/9 # How many batches have already been calculated
p <- batch*rounds

# Load null distributions
weigths_perm <- matrix(0,nrow=ncol(data_movie),ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("Movie clips: Loading null dist: ",i,"/",rounds,sep = ""))
  w <- read.csv(paste(path,"/pcoa_components_nulldist_",as.integer(i*batch),".csv",sep =""))
  weigths_perm[,((i-1)*batch+1):(i*batch)] <- as.matrix(w[,2:(batch+1)])
}

# Significance of components
s_component_clip <- c()
p_component_clip <- c()
for(j in seq(from=1, to=nrow(weigths_perm))){
  s_component_clip[j] <- sum((weigths_perm[j,] > weights_clip[j]),na.rm = TRUE)
  p_component_clip[j] <- s_component_clip[j]/p # This gives the exact p-value for each PC
}

##----------------------------------------------------------------------------------------------------------------------------------------------------
# How many significant components based on the 76 common features in full movie dataset?

path <- "/path/permutations_generalization_movie"
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

## 7 PCs were significant for movie clip data, 10 PCs for full movie data (this information is used in plotting, see figures script)

