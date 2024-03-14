# Consensus clustering analyses for studies 2, 3 & 4

# Severi Santavirta 26.2.2024

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
tbl1 <- read.csv('path/data_avg_movieclips/average_ratings_videoset_1.csv')
tbl2 <- read.csv('path/data_avg_movieclips/average_ratings_videoset_2.csv')
tbl3 <- read.csv('path/data_avg_movieclips/average_ratings_videoset_3.csv')
tbl4 <- read.csv('path/data_avg_movieclips/average_ratings_videoset_4.csv')
tbl5 <- read.csv('path/data_avg_movieclips/average_ratings_videoset_5.csv')
tbl6 <- read.csv('path/data_avg_movieclips/average_ratings_videoset_6.csv')

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
## Consensus clustering analysis of Study 2 (movie clip dataset, HC, pearson, nk=5:40, rep=1000)

# Run consensus clustering on the catenated data
data <- rbind(tbl1,tbl2,tbl3,tbl4,tbl5,tbl6)

# Run consensus clustering in reasonable amount of clusters
CC <- consensus_cluster(t(data), nk = 5:40, p.item = 0.8, reps = 1000, algorithms = c("hc"),distance = c("pearson"),scale = FALSE)

# Consensus matrix
m_trim <- consensus_matrix(CC)
rownames(m_trim) <- colnames(tbl1)
colnames(m_trim) <- colnames(tbl1)

# Save results
#save(list = c("data", "m_trim","tbl","CC"), file = "path/analysis1.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Consensus clustering analysis in Study 3, (full movie dataset, features that are common in both sets)

tbl <- read.csv("path/data_avg_fullmovie/fullmovie_ratings_common_w_movieclips.csv.csv")
data <- tbl[,2:77] 

# Run consensus clustering with all different algorithms and in reasonable amount of clusters
CC <- consensus_cluster(t(data), nk = 5:40, p.item = 0.8, reps = 1000, algorithms = c("hc"),distance = c("pearson"),scale = FALSE)

# Consensus matrix
m_trim <- consensus_matrix(CC)

names <- str_replace_all(colnames(data),"[.]"," ")
names[43] <- "Exerting self-control"
rownames(m_trim) <- names
colnames(m_trim) <- names
colnames(data) <- names

# Save results
#save(list = c("data", "m_trim","tbl","CC"), file = "path/analysis2.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Consensus clustering analysis in Study 3, (movie clip dataset, features that are common in both sets)

tbl <- read.csv("path/data_avg_fullmovie/movieclips_ratings_common_w_fullmovie.csv")
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
#save(list = c("data", "m_trim","tbl","CC"), file = "path/analysis3.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mantel test for the similarity of the clustering results for movie clip data and full movie data (Study 3)

library(ape)
library(lessR)
library(corrplot)

load("path/analysis2.RData")
data_movie <- data
m_movie <- m_trim

load("/path/analysis3.RData")
data_clips <- data
m_clips <- m_trim

# Order both matrices by the hclust order of megaperception matrix
m_clips <- corReorder(m_clips, order = "hclust", hclust_type = "average")
m_movie <- m_movie[colnames(m_clips),colnames(m_clips)]

# Test the correlation between correlation matrices
cormat_movie <- cor(data_movie)
cormat_clips <- cor(data_clips)

cormat_movie <- cormat_movie[colnames(m_clips),colnames(m_clips)]
cormat_clips <- cormat_clips[colnames(m_clips),colnames(m_clips)]

#Select the lower triangle for correlation calculation
lower_idx_movie <- lower.tri(cormat_movie)
lower_idx_clips <- lower.tri(cormat_clips)
lower_triangle_movie <- cormat_movie[lower_idx_movie]
lower_triangle_clips <- cormat_clips[lower_idx_clips]

# Mantel test, correlation matrices
r_correlation <-  cor(as.vector(lower_triangle_movie), as.vector(lower_triangle_clips))
test_cormat <- mantel.test(cormat_clips,cormat_movie,nperm = 1000000, graph = T,alternative = "greater")

#Select the lower triangle for correlation calculation
lower_idx_movie <- lower.tri(m_movie)
lower_idx_clips <- lower.tri(m_clips)
lower_triangle_movie <- m_movie[lower_idx_movie]
lower_triangle_clips <- m_clips[lower_idx_clips]

# Test the correlation between consensus matrices
r_consensus <-  cor(as.vector(lower_triangle_movie ), as.vector(lower_triangle_clips))
test_consensus <- mantel.test(m_clips,m_movie,nperm = 1000000, graph = T,alternative = "greater")

# Save results
#save(list = c("mantel_consensus_matrices", "mantel_correlation_matrices"), file = "path/mantel_test_movieclips_fullmovie.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load movie frame data of Study 4

## Read data
tbl1 <- read.csv('path/data_avg_movieframes/average_ratings_videoset_1.csv')
tbl2 <- read.csv('path/data_avg_movieframes/average_ratings_videoset_2.csv')
tbl3 <- read.csv('path/data_avg_movieframes/average_ratings_videoset_3.csv')
tbl4 <- read.csv('path/data_avg_movieframes/average_ratings_videoset_4.csv')
tbl5 <- read.csv('path/data_avg_movieframes/average_ratings_videoset_5.csv')
tbl6 <- read.csv('path/data_avg_movieframes/average_ratings_videoset_6.csv')

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
# Consensus clustering analysis of Study 4 (movie frame dataset, HC, pearson, nk=5:40, rep=1000)

# Run consensus clustering on the catenated data
data <- rbind(tbl1,tbl2,tbl3,tbl4,tbl5,tbl6)

# Run consensus clustering in reasonable amount of clusters
CC <- consensus_cluster(t(data), nk = 5:40, p.item = 0.8, reps = 1000, algorithms = c("hc"),distance = c("pearson"),scale = FALSE)

# Consensus matrix
m_trim <- consensus_matrix(CC)
rownames(m_trim) <- colnames(tbl1)
colnames(m_trim) <- colnames(tbl1)

# Save results
save(list = c("data", "m_trim","CC"), file = "/path/analysis1_movieframe.RData")

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mantel test for the similarity of the clustering results for movie clip data and movie frame data (Study 4)

library(ape)
library(lessR)
library(corrplot)

load("/path/analysis1.RData")
data_video <- data
m_video <- m_trim

load("/path/analysis1_movieframe.RData")
data_frame <- data
m_frame <- m_trim

# Order both matrices by the hclust order of megaperception matrix (drops coughing & vomiting/urinating/defecating out, since they were not perceived in the video dataset)
m_video <- corReorder(m_video, order = "hclust", hclust_type = "average")
m_frame <- m_frame[colnames(m_video),colnames(m_video)]

# Test the correlation between correlation matrices
cormat_video <- cor(data_video)
cormat_frame <- cor(data_frame)

cormat_video <- cormat_video[colnames(m_video),colnames(m_video)]
cormat_frame <- cormat_frame[colnames(m_video),colnames(m_video)]

#Select the lower triangle for correlation calculation
lower_idx_video <- lower.tri(cormat_video)
lower_idx_frame <- lower.tri(cormat_frame)
lower_triangle_video <- cormat_video[lower_idx_video]
lower_triangle_frame <- cormat_frame[lower_idx_frame]

# Mantel test, correlation matrices
r_correlation <-  cor(as.vector(lower_triangle_video), as.vector(lower_triangle_frame))
test_cormat <- mantel.test(cormat_video,cormat_frame,nperm = 1000000, graph = T,alternative = "greater")

#Select the lower triangle for correlation calculation
lower_idx_video <- lower.tri(m_video)
lower_idx_frame <- lower.tri(m_frame)
lower_triangle_video <- m_video[lower_idx_video]
lower_triangle_frame <- m_frame[lower_idx_frame]

# Test the correlation between consensus matrices
r_consensus <-  cor(as.vector(lower_triangle_video), as.vector(lower_triangle_frame))
test_consensus <- mantel.test(m_video,m_frame,nperm = 1000000, graph = T,alternative = "greater")

# Save results
#save(list = c("test_cormat", "test_consensus"), file = "/path/mantel_test_movieclips_movieframes.RData")

