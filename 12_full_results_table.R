# Collect all featurewise values into one table for the primary movie clip dataset
#       1. Feature name
#       2. Bimodality BC
#       3. Bimodality Hartigan
#       4. Occurrence rate
#       5. Standard error of the mean
#       6. Single measures ICC
#       7. Cronbach's alpha
#       7. PC component loadings 1-8 & statistical signifigance 
#       8. Cluster membership
#
# Severi Santavirta 26.2.2024

#------------------------------------------------------------------------------------------------------------------------------
library(stringr)
library(psych)
library(lessR)
library(mousetrap)

# Bimodality

# Catenate data from different videosets
data <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_1.csv")
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_2.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_3.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_4.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_5.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_6.csv"))

features <- colnames(data)
features <- str_replace_all(features,"_"," ")
features[92] <- "Exerting self-control"
tbl <- data.frame(features)

# Binomial coefficient and Hartigan's dip test (simulated P-value)
res <- mt_check_bimodality(data,methods=c("BC", "HDS_sim"))
hds_sim <- t(res$HDS_simulated_p_value)
bc <- t(res$BC)
tbl$bimodality_coefficient <- bc
tbl$hartigans_dip_test_pval <- hds_sim

# Occurrence rate, SEM, ICC:s, Cronbach's alpha
consistency_data <- read.csv("/path/consistency_measures.csv")
tbl$occurrence_rate <- consistency_data$Occurrence/234
tbl$standard_error <- consistency_data$SE
tbl$icc <- consistency_data$ICC
tbl$icc_ci95_lb <- consistency_data$ICC_lb
tbl$icc_ci95_ub <- consistency_data$ICC_ub
tbl$alpha <- consistency_data$Alpha

# PCoA feature loadings and statistical significance
# Load data
load("/path/analysis1.RData")

# Distance measure for PCoA
distance <- 1-cor(data) #Pearson

# PCoA
fit <- cmdscale(distance,eig=TRUE,k=135)
loadings <- fit$points
sign  <- read.csv("/path/pcoa_loadings_pval.csv")

tbl$pcoa_pc1_loadings <- loadings[,1]
tbl$pcoa_pc1_loadings_pval <- sign[,2]
tbl$pcoa_pc2_loadings <- loadings[,2]
tbl$pcoa_pc2_loadings_pval <- sign[,3]
tbl$pcoa_pc3_loadings <- loadings[,3]
tbl$pcoa_pc3_loadings_pval <- sign[,4]
tbl$pcoa_pc4_loadings <- loadings[,4]
tbl$pcoa_pc4_loadings_pval <- sign[,5]
tbl$pcoa_pc5_loadings <- loadings[,5]
tbl$pcoa_pc5_loadings_pval <- sign[,6]
tbl$pcoa_pc6_loadings <- loadings[,6]
tbl$pcoa_pc6_loadings_pval <- sign[,7]
tbl$pcoa_pc7_loadings <- loadings[,7]
tbl$pcoa_pc7_loadings_pval <- sign[,8]
tbl$pcoa_pc8_loadings <- loadings[,8]
tbl$pcoa_pc8_loadings_pval <- sign[,9]

# Clusters
# Load hierarchical clustering data
load("/path/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

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

clusters <- c("Introversion","Feeding","Self-control","Achievement","Emotional affection","Prosociality & Pleasant feelings & prosociality","Extraversion & playfulness","Masculinity","Social engagement","Physical discomfort","Gesturing","Emotional expression","Motivation","Antisocial behaviour","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical affection")
names <- colnames(data)
df <- data.frame(names)
df$group <- factor(cc,levels = unique(cc))

# Cluster names
cluster_names <- rep("",136)
idx <- (df$group==1 | df$group==2 | df$group==4 | df$group==5 | df$group==7 | df$group==12 | df$group==15 | df$group==17 | df$group==24 | df$group==29)
cluster_names[idx] <- "Unclustered"
cluster_names[df$group==3] <- clusters[1]
cluster_names[df$group==6] <- clusters[2]
cluster_names[df$group==8] <- clusters[3]
cluster_names[df$group==9] <- clusters[4]
cluster_names[df$group==10] <- clusters[5]
cluster_names[df$group==11] <- clusters[6]
cluster_names[df$group==13] <- clusters[7]
cluster_names[df$group==14] <- clusters[8]
cluster_names[df$group==16] <- clusters[9]
cluster_names[df$group==18] <- clusters[10]
cluster_names[df$group==19] <- clusters[11]
cluster_names[df$group==20] <- clusters[12]
cluster_names[df$group==21] <- clusters[13]
cluster_names[df$group==22] <- clusters[14]
cluster_names[df$group==23] <- clusters[15]
cluster_names[df$group==25] <- clusters[16]
cluster_names[df$group==26] <- clusters[17]
cluster_names[df$group==27] <- clusters[18]
cluster_names[df$group==28] <- clusters[19]
cluster_names[df$group==30] <- clusters[20]
df$group_name <- cluster_names

ord <- match(df$names,tbl$features)
df <- df[order(ord),]
tbl$hierarchical_cluster <- df$group_name

write.csv(tbl,file = "/path/full_result_table.csv")
