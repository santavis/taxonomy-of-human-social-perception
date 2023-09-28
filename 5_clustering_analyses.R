# Correlation analyses

library(stringr)
library(diceR)

pearson <- function(x) {
  x %>%
    t() %>%
    stats::cor(method = "pearson") %>%
    magrittr::subtract(1, .) %>%
    magrittr::extract(lower.tri(.)) %>%
    `attributes<-`(
      list(
        Size = nrow(x),
        Labels = rownames(x),
        Diag = FALSE,
        Upper = FALSE,
        methods = "pearson",
        class = "dist"
      )
    )
}

## Read data
tbl1 <- read.csv('/path/data/data_perceptual_avg/data_perceptual_avg_videoset_1.csv')
tbl2 <- read.csv('/path/data/data_perceptual_avg/data_perceptual_avg_videoset_2.csv')
tbl3 <- read.csv('/path/data/data_perceptual_avg/data_perceptual_avg_videoset_3.csv')
tbl4 <- read.csv('/path/data/data_perceptual_avg/data_perceptual_avg_videoset_4.csv')
tbl5 <- read.csv('/path/data/data_perceptual_avg/data_perceptual_avg_videoset_5.csv')
tbl6 <- read.csv('/path/data/data_perceptual_avg/data_perceptual_avg_videoset_6.csv')

# Correct dimension names
dims <- str_replace_all(colnames(tbl1),"_"," ")
dims <- str_replace_all(dims,"[.]","-")

colnames(tbl1) <- dims
colnames(tbl2) <- dims
colnames(tbl3) <- dims
colnames(tbl4) <- dims
colnames(tbl5) <- dims
colnames(tbl6) <- dims

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Analysis 1, primary dataset

# Run consensus clustering on the catenated data
data <- rbind(tbl1,tbl2,tbl3,tbl4,tbl5,tbl6)

# Run consensus clustering in reasonable amount of clusters
CC <- consensus_cluster(t(data), nk = 5:40, p.item = 0.8, reps = 1000, algorithms = c("hc"),distance = c("pearson"),scale = FALSE)

# Consensus matrix
m_trim <- consensus_matrix(CC)
rownames(m_trim) <- colnames(tbl1)
colnames(m_trim) <- colnames(tbl1)

# Save results
#save(list = c("data", "m_trim","tbl","CC"), file = "/path/consensus_clustering_analysis/analyses/analysis1.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Analysis 2, validation dataset

tbl <- read.csv("/path/data/data_perceptual_avg/validation_set/ratings_common_w_primary_set.csv")
data <- tbl[,2:77] 

# Run consensus clustering in reasonable amount of clusters
CC <- consensus_cluster(t(data), nk = 5:40, p.item = 0.8, reps = 1000, algorithms = c("hc"),distance = c("pearson"),scale = FALSE)

# Consensus matrix
m_trim <- consensus_matrix(CC)

names <- str_replace_all(colnames(data),"[.]"," ")
names[43] <- "Exerting self-control"
rownames(m_trim) <- names
colnames(m_trim) <- names
colnames(data) <- names

# Save results
#save(list = c("data", "m_trim","tbl","CC"), file = "/path/consensus_clustering_analysis/analyses/analysis2.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Analysis 3, primary dataset, features common with validation dataset

tbl <- read.csv("/path/data/data_perceptual_avg/ratings_common_w_validation_set.csv")
data <- tbl[,2:77] 

# Run consensus clustering with all different algorithms and in reasonable amount of clusters
CC <- consensus_cluster(t(data), nk = 5:40, p.item = 0.8, reps = 1000, algorithms = c("hc"),distance = c("pearson"),scale = FALSE)

# Consensus matrix
m_trim <- consensus_matrix(CC)

names <- str_replace_all(colnames(data),"[.]"," ")
names[50] <- "Exerting self-control"
rownames(m_trim) <- names
colnames(m_trim) <- names
colnames(data) <- names

# Save results
#save(list = c("data", "m_trim","tbl","CC"), file = "/path/consensus_clustering_analysis/analyses/analysis3.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mantel test for matrices (analysis 2 and analysis 3)

library(ape)
library(lessR)
library(corrplot)

load("/path/consensus_clustering_analysis/analyses/analysis2.RData")
data_movie <- data
m_movie <- m_trim

load("/path/consensus_clustering_analysis/analyses/analysis3.RData")
data_clips <- data
m_clips <- m_trim

# Order both matrices by the hclust order of primary dataset matrix
m_clips <- corReorder(m_clips, order = "hclust", hclust_type = "average")
m_movie <- m_movie[colnames(m_clips),colnames(m_clips)]

# Test the correlation between correlation matrices
cormat_movie <- cor(data_movie)
cormat_clips <- cor(data_clips)

cormat_movie <- cormat_movie[colnames(m_clips),colnames(m_clips)]
cormat_clips <- cormat_clips[colnames(m_clips),colnames(m_clips)]

# Mantel test, correlation matrices
r_correlation <-  cor(as.vector(cormat_clips), as.vector(cormat_movie))
test_cormat <- mantel.test(cormat_clips,cormat_movie,nperm = 1000000, graph = T,alternative = "greater")

# Test the correlation between consensus matrices
r_consensus <-  cor(as.vector(m_clips), as.vector(m_movie))
test_consensus <- mantel.test(m_clips,m_movie,nperm = 1000000, graph = T,alternative = "greater")

# Save results
#save(list = c("mantel_consensus_matrices", "mantel_correlation_matrices"), file = "/path/consensus_clustering_analysis/analyses/mantel_test_generalization.RData")
