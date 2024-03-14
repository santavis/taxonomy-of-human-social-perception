# Analyze individual feature results (ICCs, and bimodality coefficients) between the movie clip dataset and the static movie frame dataset (Study 4)

# Severi Santavirta 26.2.2024

library(ggplot2)
library(mousetrap)
library(stringr)

data_clips <- read.csv("/scratch/severi/projects/megaperception/consistency_measures.csv")
data_frames <- read.csv("/scratch/severi/projects/megaperception/consistency_measures_movieframe.csv")

# Investigate ICC
icc_mu_clips <- mean(data_clips$ICC)
icc_min_clips <- min(data_clips$ICC)
icc_max_clips <- max(data_clips$ICC)
icc_mu_frames <- mean(data_frames$ICC)
icc_min_frames <- min(data_frames$ICC)
icc_max_frames <- max(data_frames$ICC)
icc_cor <- cor.test(data_clips$ICC,data_frames$ICC)
icc_diff <- data_clips$ICC - data_frames$ICC
icc_ttest <- t.test(icc_diff)
hist(icc_diff)

# Investigate bimodality
# Catenate data from different videosets (dynamic dataset)
data_clips <- read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieclips/average_ratings_videoset_1.csv")
data_clips <- rbind(data_clips,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieclips/average_ratings_videoset_2.csv"))
data_clips <- rbind(data_clips,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieclips/average_ratings_videoset_3.csv"))
data_clips <- rbind(data_clips,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieclips/average_ratings_videoset_4.csv"))
data_clips <- rbind(data_clips,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieclips/average_ratings_videoset_5.csv"))
data_clips <- rbind(data_clips,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieclips/average_ratings_videoset_6.csv"))

# Catenate data from different videosets (static dataset)
data_frames <- read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_1.csv")
data_frames <- rbind(data_frames,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_2.csv"))
data_frames <- rbind(data_frames,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_3.csv"))
data_frames <- rbind(data_frames,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_4.csv"))
data_frames <- rbind(data_frames,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_5.csv"))
data_frames <- rbind(data_frames,read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_6.csv"))

# Binomial coefficient and Hartigan's dip test
res_clips <- mt_check_bimodality(data_clips,methods=c("BC", "HDS_sim"))
res_frames <- mt_check_bimodality(data_frames,methods=c("BC", "HDS_sim"))

# Compare BC
bc_clips <- as.data.frame(t(res_clips$BC))
colnames(bc_clips) <- "bc"
bc_clips$feature <- rownames(bc_clips)
bc_frames <- as.data.frame(t(res_frames$BC))
colnames(bc_frames) <- "bc"
bc_frames$feature <- rownames(bc_frames) 
bc_corr <- cor.test(bc_clips$bc,bc_frames$bc)
bc_sign_clips <- bc_clips[bc_clips$bc>0.555,] 
bc_sign_frames <- bc_frames[bc_frames$bc>0.555,]
bc_intersect <- intersect(bc_sign_clips$feature,bc_sign_frames$feature)

# Compare HDS
hds_clips <- as.data.frame(t(res_clips$HDS_simulated_p_value))
colnames(hds_clips) <- "hds"
hds_clips$feature <- rownames(hds_clips)
hds_frames <- as.data.frame(t(res_frames$HDS_simulated_p_value))
colnames(hds_frames) <- "hds"
hds_frames$feature <- rownames(hds_frames) 
hds_corr <- cor.test(hds_clips$hds,hds_frames$hds)
hds_sign_clips <- hds_clips[bc_clips$bc>0.555 & hds_clips$hds<0.05,] 
hds_sign_frames <- hds_frames[bc_frames$bc>0.555 & hds_frames$hds<0.05,]
hds_intersect <- intersect(hds_sign_clips$feature,hds_sign_frames$feature)


