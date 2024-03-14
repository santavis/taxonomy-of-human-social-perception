# Calculate how many significant components movie clip data (136 features) and movie frame data (136 features) have. (Study 4)
# This information is used when comparing the component loadings across datasets (see figures script)

# Severi Santavirta 26.2.2024

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate the null distributions for component weights for movie frame data (primary movie clip component significances have been calculated previously)

library(stringr)

# Load frame data  
load("/path/analysis1_movieframe.RData")
data_frame <- data

# Distance measure
distance_frame <- 1-cor(data_frame) #Pearson

# PCoA with real data
fit_frame <- cmdscale(distance_frame,eig=TRUE,k=135)
loadings_frame <- fit_frame$points
dim_names <- colnames(data_frame)

# Run permutations on randomized data to estimate chance distributions

p <- 1000000 # How many permutations
batch <- 1000
n <- 0 # set own seed for each sample
k <- 0 # Count permutations
weights_perm <- matrix(0,nrow = ncol(data_frame), ncol = batch)

for(i in seq(from=1,to=p)){
  k <- k+1
  if(i%%100 == 0){
    print(paste("Movie: permutation: ",i,"/",p))
  }
  data_perm <- data_frame
  for(j in seq(from=1,to=ncol(data_perm))){ # Create bootsrapped sample
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
    #write.csv(weights_perm,paste("/path/permutations_movieframe/pcoa_components_nulldist_",i,".csv",sep = ""))
    k<-0
  }
}


rm(list = ls())


##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate how many significant components movie clip data and movie frame data have. (Study 4)
# This information is used when comparing the component loadings across datasets (see figures script)
#
# Severi Santavrita 26.2.2024

# PCoA generalization

library(stringr)
library(corrplot)
library(lessR)

# Calculate the real PCoA solutions first for both datasets

# Load average data from the clustering datasets
load("/path/analysis1.RData")
data_clips <- data

load("/path/analysis1_movieframe.RData")
data_frames <- data

# Distance measure
distance_clips <- 1-cor(data_clips) #Pearson
distance_frames <- 1-cor(data_frames) #Pearson

# PCoA with real data
fit_clips <- cmdscale(distance_clips,eig=TRUE,k=135)
fit_frames <- cmdscale(distance_frames,eig=TRUE,k=135)
loadings_clips <- fit_clips$points
loadings_frames <- fit_frames$points
weights_frames <- fit_frames$eig
weights_clips <- fit_clips$eig

##------------------------------------------------------------------------------
# How many significant components based on the movie clip data?
path <- "/path/permutations"
batch <- 1000 # Batch size
rounds <- length(dir(path))/9 # How many batches have already been calculated
p <- batch*rounds

# Load null distributions
weigths_perm <- matrix(0,nrow=ncol(data_clips),ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("Clips: Loading null dist: ",i,"/",rounds,sep = ""))
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

##------------------------------------------------------------------------------
# How many significant components based on static frames?
path <- "/path/permutations_movieframe"
batch <- 1000 # Batch size
rounds <- length(dir(path)) # How many batches have already been calculated
p <- batch*rounds

# Load null distributions
weigths_perm <- matrix(0,nrow=ncol(data_frames),ncol=batch*rounds)

for(i in seq(from=1,to=rounds)){
  print(paste("frames: Loading null dist: ",i,"/",rounds,sep = ""))
  w <- read.csv(paste(path,"/pcoa_components_nulldist_",as.integer(i*batch),".csv",sep =""))
  weigths_perm[,((i-1)*batch+1):(i*batch)] <- as.matrix(w[,2:(batch+1)])
}

# Significance of components
s_component_frames <- c()
p_component_frames <- c()
for(j in seq(from=1, to=nrow(weigths_perm))){
  s_component_frames[j] <- sum((weigths_perm[j,] > weights_frames[j]),na.rm = TRUE)
  p_component_frames[j] <- s_component_frames[j]/p # This gives the exact p-value for each PC
}

## 8 PCs were significant for movie clip data, 9 PCs for movie frame data (this information is used in plotting, see figures script)

