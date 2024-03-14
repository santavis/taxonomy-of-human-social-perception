# ICC stability analysis of ratings (Two way random, absolute agreement, single measures)
# This analysis is used to estimate how stable ICC estimates are with different amount of subject to choose how much data should be collected. 
# A previously collected dataset of emotion ratings is used as reference
#
#   1. Stability of single measures ICC in emotion features (emotion data not publicly available)
#   3. Stability of single measures ICC in social features 
#
# Severi Santavirta 26.2.2024
#
#---------------------------------------------------------------------------------------------------------------------------------------------
# 1. Stability of single measures ICC in emotion features 

library(dplyr)
library(stringr)
library(psych)

# How many combination to draw for ICC estimation
perm <- 1000

# Path to save data
path <- "/path/results"

# Load emotion data (not publicly available)
data_wide <- read.csv('/path/emotion_full_data_wide.csv')
cols <- colnames(data_wide)
dims <- cols[7:39] 

for(i in seq(from=7,to=length(dims))){
  
  # Create folder
  dir.create(paste(path,"/",dim_read[i],sep = ""))
  dir.create(paste(path,"/",dim_read[i],"/distributions",sep = ""))
  
  vset1 <- c()
  vset2 <- c()
  vset3 <- c()
  
  data_dim <- data_wide[,c(1,4,5,i)]
  data_dim <- data_dim[!is.na(data_dim[,4]),] # Delete empty rows
  colnames(data_dim) <- c("subject","videoset","video","rating")
  
  subjects <- unique(data_dim$subject)
  for(j in seq(from=1,to=length(subjects))){
    data_sub <- data_dim[data_dim$subject==subjects[j],]
    data_sub <- data_sub[order(data_sub$video),]
    colnames(data_sub) <- c("subject","videoset","video",subjects[j])
    
    if(data_sub$videoset[1] == "VideoZoneSubset1") {
      vset1 <- cbind(vset1,data_sub[,4])
    }
    if(data_sub$videoset[1] == "VideoZoneSubset2") {
      vset2 <- cbind(vset2,data_sub[,4])
    }
    if(data_sub$videoset[1] == "VideoZoneSubset3") {
      vset3 <- cbind(vset3,data_sub[,4])
    }
  }
  
  r_idx <- cumsum(c(nrow(vset1),nrow(vset2),nrow(vset3)))
  c_idx <- cumsum(c(ncol(vset1),ncol(vset2),ncol(vset3)))
  
  # Create matrix for ICC calculations
  m <- matrix(ncol = c_idx[3],nrow = r_idx[3])
  m[1:r_idx[1],1:c_idx[1]] <- vset1
  m[(r_idx[1]+1):r_idx[2],(c_idx[1]+1):c_idx[2]] <- vset2
  m[(r_idx[2]+1):r_idx[3],(c_idx[2]+1):c_idx[3]] <- vset3
  
  # Shortest set
  min_subs <- min(ncol(vset1),ncol(vset2),ncol(vset3))
  
  dim_icc_lb <- rep(0, min_subs-1)
  dim_icc_ub <- rep(0, min_subs-1)
  dim_ci_lb <- rep(0, min_subs-1)
  dim_ci_ub <- rep(0, min_subs-1)
  dim_perms <- rep(0, min_subs-1)
  n <- 0
  for(j in seq(from=2, to=min_subs)){ # Number of subjects in a set
    print(paste("dim:",i,"set",j))
    
    # Sample random sets of subjects from each videoset for ICC calculation
    v1_set <- combn(seq(from=1,to=c_idx[1]),j)
    v2_set <- combn(seq(from=c_idx[1]+1,to=c_idx[2]),j)
    v3_set <- combn(seq(from=c_idx[2]+1,to=c_idx[3]),j)
    
    # Maximum number of combination
    max_perm <- ncol(v1_set)*ncol(v2_set)*ncol(v3_set)
    
    if(is_na(max_perm)){ # Very large number of possible combinations
      p <- perm
    } else {
      if(max_perm > perm){ # Very large number of possible combinations
        p <- perm
      } else { # Only few possible combinations
        p <- max_perm
      }
    }
    
    # For each permutation and set within a dimension there is a different seed, but for different dimensions the seeds are the same
    # so if there is some bias the bias is the same for each social dimension
    set.seed(n+1)
    perm1 <- sample(seq(from=1,to=ncol(v1_set)),p,replace = TRUE)
    set.seed(n+2)
    perm2 <- sample(seq(from=1,to=ncol(v2_set)),p,replace = TRUE)
    set.seed(n+3)
    perm3 <- sample(seq(from=1,to=ncol(v3_set)),p,replace = TRUE)
    n <- n+3
    
    # Columns of the combination matrices
    ind_set <- cbind(perm1,perm2,perm3)
    
    # Permute ICC
    perm_lb <- c()
    perm_icc <- c()
    perm_ub <- c()
    for(l in seq(from=1, to=p)){ # permutations
      
      # Columns of the data matrix
      ind_col <- c(v1_set[,ind_set[l,1]],
                   v2_set[,ind_set[l,2]],
                   v3_set[,ind_set[l,3]])
      
      s <- m[,ind_col]
      res <- ICC(s,alpha=.05)
      perm_lb <- c(perm_lb,res$results$`lower bound`[2])
      perm_icc <- c(perm_icc,res$results$ICC[2])
      perm_ub <- c(perm_ub,res$results$`upper bound`[2])
    }
    
    # Save permuted distributions
    df <- data.frame(perm_icc,perm_lb,perm_ub)
    write.csv(df,paste(path,"/",dims[i],"/distributions/",j,"_subjects_single_measures.csv",sep = ""))
    
    # To assess how confident we are with the ICC, we ...
    #   1. calculate the 2.5% and 97.5% quantiles of permuted ICC values
    #   2. calculate the 2.5% quantile of permuted lower CI and 97.5% quantile of permuted upper CI
    dim_icc_lb[j-1] <- quantile(perm_icc,0.025)
    dim_icc_ub[j-1] <- quantile(perm_icc,0.975)
    dim_ci_lb[j-1] <- quantile(perm_lb,0.025)
    dim_ci_ub[j-1] <- quantile(perm_ub,0.975)
    dim_perms[j-1] <- p
    
  }
  df <- data.frame(dim_icc_lb,dim_icc_ub,dim_ci_lb,dim_ci_ub,dim_perms)
  write.csv(df,paste(path,"/",dims[i],"/results_single_measures.csv",sep = ""))
}

#---------------------------------------------------------------------------------------------------------------------------------------------
# 3. Stability of single measures ICC in social features 

library(dplyr)
library(stringr)
library(psych)

# How many combination to draw for ICC estimation
perm <- 1000

# Path to save data
path <- "/path/icc_analysis"

# Catenate data from different videosets
data <- read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_1.csv")
data <- rbind(data,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_2.csv"))
data <- rbind(data,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_3.csv"))
data <- rbind(data,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_4.csv"))
data <- rbind(data,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_5.csv"))
data <- rbind(data,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_6.csv"))

# Dim names
dims <- colnames(data)
dims <- str_replace_all(dims,"_"," ")
dims <- str_replace_all(dims,"[.]","-")

dim_read <- colnames(data)
dim_read[92] <- "Exerting_self-control"

datapath <- "/scratch/severi/projects/megaperception/data/data_git/data_csv_movieclips"

for(i in seq(from=1, to=length(data))){ # Dimensions
  
  # Create folder
  dir.create(paste(path,"/",dim_read[i],sep = ""))
  dir.create(paste(path,"/",dim_read[i],"/distributions",sep = ""))
  
  # Load data
  v1 <- read.csv(paste(datapath,"/",dim_read[i],"_1.csv",sep = ""))
  v2 <- read.csv(paste(datapath,"/",dim_read[i],"_2.csv",sep = ""))
  v3 <- read.csv(paste(datapath,"/",dim_read[i],"_3.csv",sep = ""))
  v4 <- read.csv(paste(datapath,"/",dim_read[i],"_4.csv",sep = ""))
  v5 <- read.csv(paste(datapath,"/",dim_read[i],"_5.csv",sep = ""))
  v6 <- read.csv(paste(datapath,"/",dim_read[i],"_6.csv",sep = ""))
  
  # Indices
  idx <- cumsum(c(length(v1),length(v2),length(v3),length(v4),length(v5),length(v6)))
  
  # Create a combined design matrix for ICC (absolute) calculation over all videosets
  m <- matrix(NA,234,idx[6])
  m[1:39,(1:idx[1])] <- as.matrix(v1)
  m[40:78,((idx[1]+1):idx[2])] <- as.matrix(v2)
  m[79:117,((idx[2]+1):idx[3])] <- as.matrix(v3)
  m[118:156,((idx[3]+1):idx[4])] <- as.matrix(v4)
  m[157:195,((idx[4]+1):idx[5])] <- as.matrix(v5)
  m[196:234,((idx[5]+1):idx[6])] <- as.matrix(v6)
  
  # Shortest set
  min_subs <- min(length(v1),length(v2),length(v3),length(v4),length(v5),length(v6))
  
  dim_icc_lb <- rep(0, min_subs-1)
  dim_icc_ub <- rep(0, min_subs-1)
  dim_ci_lb <- rep(0, min_subs-1)
  dim_ci_ub <- rep(0, min_subs-1)
  dim_perms <- rep(0, min_subs-1)
  n <- 0
  for(j in seq(from=2, to=min_subs)){ # Number of subjects in a set
    print(paste("dim:",i,"set",j))
    
    # Sample random sets of subjects from each videoset for ICC calculation
    v1_set <- combn(seq(from=1,to=idx[1]),j)
    v2_set <- combn(seq(from=idx[1]+1,to=idx[2]),j)
    v3_set <- combn(seq(from=idx[2]+1,to=idx[3]),j)
    v4_set <- combn(seq(from=idx[3]+1,to=idx[4]),j)
    v5_set <- combn(seq(from=idx[4]+1,to=idx[5]),j)
    v6_set <- combn(seq(from=idx[5]+1,to=idx[6]),j)
    
    # Maximum number of combination
    max_perm <- ncol(v1_set)*ncol(v2_set)*ncol(v3_set)*ncol(v4_set)*ncol(v5_set)*ncol(v6_set)
    
    if(is_na(max_perm)){ # Very large number of possible combinations
      p <- perm
    } else {
      if(max_perm > perm){ # Very large number of possible combinations
        p <- perm
      } else { # Only few possible combinations
        p <- max_perm
      }
    }
    
    # For each permutation and set within a dimension there is a different seed, but for different dimensions the seeds are the same
    # so if there is some bias the bias is the same for each social dimension
    set.seed(n+1)
    perm1 <- sample(seq(from=1,to=ncol(v1_set)),p,replace = TRUE)
    set.seed(n+2)
    perm2 <- sample(seq(from=1,to=ncol(v2_set)),p,replace = TRUE)
    set.seed(n+3)
    perm3 <- sample(seq(from=1,to=ncol(v3_set)),p,replace = TRUE)
    set.seed(n+4)
    perm4 <- sample(seq(from=1,to=ncol(v4_set)),p,replace = TRUE)
    set.seed(n+5)
    perm5 <- sample(seq(from=1,to=ncol(v5_set)),p,replace = TRUE)
    set.seed(n+6)
    perm6 <- sample(seq(from=1,to=ncol(v6_set)),p,replace = TRUE)
    n <- n+6
    
    # Columns of the combination matrices
    ind_set <- cbind(perm1,perm2,perm3,perm4,perm5,perm6)
    # Permute ICC
    perm_lb <- c()
    perm_icc <- c()
    perm_ub <- c()
    for(l in seq(from=1, to=p)){ # permutations
      
      # Columns of the data matrix
      ind_col <- c(v1_set[,ind_set[l,1]],
                   v2_set[,ind_set[l,2]],
                   v3_set[,ind_set[l,3]],
                   v4_set[,ind_set[l,4]],
                   v5_set[,ind_set[l,5]],
                   v6_set[,ind_set[l,6]])
      
      s <- m[,ind_col]
      res <- ICC(s,alpha=.05)
      perm_lb <- c(perm_lb,res$results$`lower bound`[2])
      perm_icc <- c(perm_icc,res$results$ICC[2])
      perm_ub <- c(perm_ub,res$results$`upper bound`[2])
    }
    
    # Save permuted distributions
    df <- data.frame(perm_icc,perm_lb,perm_ub)
    write.csv(df,paste(path,"/",dim_read[i],"/distributions/",j,"_subjects_single_measures.csv",sep = ""))
    
    # To assess how confident we are with the ICC, we ...
    #   1. calculate the 2.5% and 97.5% quantiles of permuted ICC values
    #   2. calculate the 2.5% quantile of permuted lower CI and 97.5% quantile of permuted upper CI
    dim_icc_lb[j-1] <- quantile(perm_icc,0.025)
    dim_icc_ub[j-1] <- quantile(perm_icc,0.975)
    dim_ci_lb[j-1] <- quantile(perm_lb,0.025)
    dim_ci_ub[j-1] <- quantile(perm_ub,0.975)
    dim_perms[j-1] <- p
    
  }
  df <- data.frame(dim_icc_lb,dim_icc_ub,dim_ci_lb,dim_ci_ub,dim_perms)
  write.csv(df,paste(path,"/",dim_read[i],"/results_single_measures.csv",sep = ""))
}
