# Principal coordinate analysis (PCoA) of social perceptual features and null distributions for the permutation testing in the primary movie clip dataset (Study 2)

# Severi Santavirta 26.2.2024

library(stringr)

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

# Run permutations on randomized data to estimate chance distributions

p <- 1000000 # How many permutations
batch <- 1000
n <- 0 # set own seed for each sample
k <- 0 # Count permutations
weights_perm <- matrix(0,nrow = 136, ncol = batch)
loadings_perm1 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm2 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm3 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm4 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm5 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm6 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm7 <- matrix(0,nrow = 136, ncol = batch)
loadings_perm8 <- matrix(0,nrow = 136, ncol = batch)

for(i in seq(from=1,to=p)){
  k <- k+1
  if(i%%100 == 0){
    print(paste("Permutation: ",i,"/",p))
  }
  data_perm <- data
  for(j in seq(from=1,to=136)){ # Create boostrapped sample
    n <- n+1
    set.seed(n)
    data_perm[,j] <- sample(data_perm[,j],replace = FALSE) # Selects all the data in a random order (based on the seed number, which is different for each column and each permutation)
  }
  distance_perm <- 1-cor(data_perm)
  fit_perm <- cmdscale(distance_perm,eig=TRUE,k=135)
  
  # We know that only first 8 components are significant so store the loadings only for these components
  loadings_perm <- fit_perm$points
  loadings_perm1[,k] <- loadings_perm[,1]
  loadings_perm2[,k] <- loadings_perm[,2]
  loadings_perm3[,k] <- loadings_perm[,3]
  loadings_perm4[,k] <- loadings_perm[,4]
  loadings_perm5[,k] <- loadings_perm[,5]
  loadings_perm6[,k] <- loadings_perm[,6]
  loadings_perm7[,k] <- loadings_perm[,7]
  loadings_perm8[,k] <- loadings_perm[,8]
  
  # Store weights to test the components significance
  weights_perm[,k] <- fit_perm$eig
  
  if(i%%batch == 0){ # Save null distributions for every batch
    
    rownames(weights_perm) <- dim_names
    rownames(loadings_perm1) <- dim_names
    rownames(loadings_perm2) <- dim_names
    rownames(loadings_perm3) <- dim_names
    rownames(loadings_perm4) <- dim_names
    rownames(loadings_perm5) <- dim_names
    rownames(loadings_perm6) <- dim_names
    rownames(loadings_perm7) <- dim_names
    rownames(loadings_perm8) <- dim_names
    
    #write.csv(weights_perm,paste("/path/permutations/pcoa_components_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm1,paste("/path/permutations/pcoa_loadings_pc1_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm2,paste("/path/permutations/pcoa_loadings_pc2_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm3,paste("/path/permutations/pcoa_loadings_pc3_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm4,paste("/path/permutations/pcoa_loadings_pc4_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm5,paste("/path/permutations/pcoa_loadings_pc5_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm6,paste("/path/permutations/pcoa_loadings_pc6_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm7,paste("/path/permutations/pcoa_loadings_pc7_nulldist_",i,".csv",sep = ""))
    #write.csv(loadings_perm8,paste("/path/permutations/pcoa_loadings_pc8_nulldist_",i,".csv",sep = ""))
    
    k<-0
  }
}