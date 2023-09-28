# Supplementary Figures / materials
#
#       1. Bimodality analysis and plot (FigSI1)
#       2. Heatmap of social feature ratings (FigSI2)
#       3. Concordance analysis, polar plots of main HC dimensions (FigSI3)
#       4. Representative clustering results of the two datasets (FigSI4)
#       5. Estimate sufficient number of raters by ICC stability (FigSI5)
#       6. 3D htmls of the PCoA loadings (Supplementary media)


##-----------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Calculate and plot bimodality for the average ratings (FigSI1)

library(ggplot2)
library(mousetrap)
library(stringr)

# Catenate data from different videosets
data <- read.csv("/path/data/data_perceptual_avg/data_perceptual_avg_videoset_1.csv")
data <- rbind(data,read.csv("/path/data/data_perceptual_avg/data_perceptual_avg_videoset_2.csv"))
data <- rbind(data,read.csv("/path/data/data_perceptual_avg/data_perceptual_avg_videoset_3.csv"))
data <- rbind(data,read.csv("/path/data/data_perceptual_avg/data_perceptual_avg_videoset_4.csv"))
data <- rbind(data,read.csv("/path/data/data_perceptual_avg/data_perceptual_avg_videoset_5.csv"))
data <- rbind(data,read.csv("/path/data/data_perceptual_avg/data_perceptual_avg_videoset_6.csv"))

# Binomial coefficient and Hartigan's dip test
res <- mt_check_bimodality(data,methods=c("BC", "HDS_sim"))

# Create a data frame for plotting
stat <- data.frame(t(res$BC)[,1])
stat$hds_sim <- t(res2$HDS_simulated_p_value)[,1]
stat$name <- str_replace_all(colnames(data),"_"," ")
colnames(stat) <- c("bc","hds_sim","names")

# Create significance asterisks
sign <- rep(" ",136)
sign[stat$bc>0.555 & stat$hds_sim<0.05] <- "**"
sign[stat$bc>0.555 & !stat$hds_sim<0.05] <- "*"

# Barplot sorted Bimodality coeffs witth significance 
#png(file = "/path/visualizations/FigSI1_bimodality/FigSI1_bmodality_bars_bc.png",height = 7500,width = 5000)
ggplot(stat,aes(x = reorder(names, bc), y = bc)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x= reorder(names, bc),y=bc,label = sign), position = position_dodge(width = 1),vjust = 0.8, hjust = 0, size = 60,inherit.aes = TRUE) +
  theme_minimal() +
  ggtitle("Bimodality coefficient (BC)") +
  theme(plot.title = element_text(hjust = 0.5,size = 140),
        axis.text.y = element_text(vjust = 0.4, hjust=1, size = 70),
        axis.text.x = element_text(size = 70),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(0,0.95), expand = c(0,0))
dev.off()

##-----------------------------------------------------------------------------------------------------------------------------------------------------
## 2. Superheat population average ratings to visualize the variation in the ratings (FigSI2)

library(superheat)
library(lessR)

# Load and order data
load("/path/analyses/analysis1a.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

#pdf("/path/visualizations/FigSI2_superheat_ratings/superheat.pdf",height = 14,width = 16)
par(mar=c(2,10,2,2))
superheat(t(data),
          grid.hline = FALSE,
          force.left.label = TRUE,
          left.label.size = 0.3,
          left.label.text.size = 2.5,
          left.label.text.alignment = "right",
          left.label.col = "white",
          heat.pal = c("white","#F4A582","#B2182B"),
          heat.lim = c(0,100),
          legend = T)
#dev.off()

##---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Concordance analysis: Polarplot HC clusters based on their PCoA loadings (FigSI3)

library(ggplot2)
library(lessR)
library(corrplot)

# Load hierarchical clustering data and correct names
load("/path/analyses/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]
dims <- colnames(data)
dims <- str_replace_all(dims,"_"," ")

# Distance measure for PCoA
distance <- 1-cor(data) #Pearson

# PCoA
fit <- cmdscale(distance,eig=TRUE,k=135) # Only 56 first components have positive eigenvalues
loadings <- fit$points

# Cluster membership by cluster (Hard-coded cluster indices from consensus hierarchical clustering in Figure 54)
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

# Data frame for plotting
names <- colnames(data)
df <- data.frame(names,loadings)
df$group <- factor(cc,levels = unique(cc))

# Exclude unclustered features from the plot 
idx <- (df$group==1 | df$group==2 | df$group==4 | df$group==5 | df$group==7 | df$group==12 | df$group==15 | df$group==17 | df$group==24 | df$group==29)
df1 <- df[!idx,]
names1 <- names[!idx]
idx <- unique(df1$group)

# Load cluster significance information
p_cluster <- read.csv("/path/pcoa/pcoa_cluster_loadings_pval.csv")
sign001_cluster <- read.csv("/path/pcoa/pcoa_cluster_loadings_pval_sing001.csv")
signbonf_cluster <- read.csv("/path/pcoa/pcoa_cluster_loadings_pval_signbonf.csv")
colnames(p_cluster) <-c("PC","Introversion","Feeding","Excecutive control","Comptence & ambition","Affection","Prosociality & pleasant feelings","Extraversion & humor","Masculinity","Verbal communication","Sympathetic activity","Gesturing","Nonverbal communication","Goal oriented behaviour","Antisociality","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical interaction")

#Create significance asterisk table
sign <- character(length = 8*20)
k<-0
for(i in seq(from=1,to=20)){
  for(j in seq(from=1,to=8)){
    k<-k+1
    if(signbonf_cluster[j,(i+1)]){
      sign[k] <- "**"
    } else {
      if(sign001_cluster[j,(i+1)]){
        sign[k] <- "*"
      }
    }
  }
}
sign <- matrix(sign,nrow=8,ncol=20,byrow=FALSE)

# Main dimension labels for plot
pcoa <- as.factor(c("Pleasant - Unpleasant","Submission - Dominance","Impulsive - Cognitive","Careless - Conscientious","Introversion - Extraversion","Friendship - Sexuality","Alone - With others","Femininity - Masculinity"))
pcoa_neg <- c("Pleasant","Submission","Cognitive","Conscientious","Extraversion","Friendship","Alone","Masculinity")
pcoa_pos <- c("Unpleasant","Dominance","Impulsive","Careless","Introversion","Sexuality","With others","Femininity")
clusters <- c("Introversion","Feeding","Excecutive control","Comptence & ambition","Affection","Prosociality & pleasant feelings","Extraversion & humor","Masculinity","Verbal communication","Sympathetic activity","Gesturing","Nonverbal communication","Goal oriented behaviour","Antisociality","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical interaction")
d <- as.data.frame(pcoa)

# Polarplot main dimension loadings for each HC cluster
cluster_loadings <- c()
require(lattice)
for(i in seq(1:20)){
  d$loadings <- colMeans(as.matrix(df1[df1$group==idx[i],2:9]))
  cluster_loadings <- rbind(cluster_loadings,d$loadings)
  d$signi <- sign[,i] 
  
  # Polar plot (direction colorcoded)
  #pdf(file = paste("/path/visualizations/FigSI3_polarplots/polars_noname/",clusters[i],".pdf",sep = ""),width = 10,height=10)
  d$loadings_abs <- abs(d$loadings)
  p <- ggplot(d,aes(x=pcoa,y=loadings_abs)) +
    geom_segment(aes(x = pcoa, y = 0, xend = pcoa, yend = loadings_abs),size = abs(d$loadings)*10,arrow = arrow(length = unit(0.5, "cm")), color = ifelse((d$loadings>0),"#D6604D", "#4393C3"),stat = "identity") +
    scale_x_discrete(limits=d$pcoa) +
    geom_text(aes(x=pcoa,y=loadings_abs+0.05,label = signi),size = 14) +
    ylim(0,0.75) +
    coord_polar() +
    ggtitle(clusters[i]) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_line(size = 1),
          axis.title.x = element_blank())
  print(p)
  #dev.off()
}

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. Representative clustering of the two datasets (movie clips / full movie) that serves as the ground for the alluvial diagram in Figure 6 (FigSI4)

library(corrplot)

load("/path/clustering_analyses/analysis2.RData")
m_movie <- m_trim

load("/path/clustering_analyses/analysis3.RData")
m_clips <- m_trim

# Choose k for clips
#pdf("/path/visualizations/FigSI4_generalization_consensusmatrices/materials/consensus_clips_k16.pdf",width = 22,height = 20)
corrplot(m_clips, addrect = 16, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),)
#dev.off()

# Choose k for movie
#pdf("/path/visualizations/FigSI4_generalization_consensusmatrices/materials/consensus_movie_k19.pdf",width = 22,height = 20)
corrplot(m_movie, addrect = 19, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
#dev.off()

##-----------------------------------------------------------------------------------------------------------------------------------------------------
## 5. ICC stability plot (emotion data is not publicly shared) (Figure SI5)

library(ggplot2)
library(stringr)

# Emotion dimensions
emotion_dims <- dir("/path/icc_analysis_emotions/results") # Path to appraisals data
emotion_dim_names <- str_replace_all(emotion_dims,"_"," ")

# Social dimensions
social_dims <- dir("/path/icc_analysis/results")
social_dim_names <- str_replace_all(social_dims,"_"," ")


# Read ICC data for emotions
for(i in seq(from=1,to=length(emotion_dims))){
  
  # Emotion single measures ICC
  tbl <- read.csv(paste("/path/icc_analysis_emotions/results/",emotion_dims[i],"/results_single_measures.csv",sep=""))
  tbl$dim <- rep(emotion_dim_names[i],nrow(tbl))
  if(i==1){
    df_single_emotion <- tbl
  } else {
    df_single_emotion <- rbind(df_single_emotion,tbl)
  }
  
  # Emotion average measures ICC
  tbl <- read.csv(paste("/path/icc_analysis_emotions/results/",emotion_dims[i],"/results_average_measures.csv",sep=""))
  tbl$dim <- rep(emotion_dim_names[i],nrow(tbl))
  if(i==1){
    df_average_emotion <- tbl
  } else {
    df_average_emotion <- rbind(df_average_emotion,tbl)
  }
}

# Read ICC data for social features
for(i in seq(from=1,to=length(social_dims))){
  
  # Social single measures ICC
  tbl <- read.csv(paste("/path/icc_analysis/results/",social_dims[i],"/results_single_measures.csv",sep=""))
  tbl$dim <- rep(social_dim_names[i],nrow(tbl))
  if(i==1){
    df_single_social <- tbl
  } else {
    df_single_social <- rbind(df_single_social,tbl)
  }
  
  # Social average measures ICC
  tbl <- read.csv(paste("/path/icc_analysis/results/",social_dims[i],"/results_average_measures.csv",sep=""))
  tbl$dim <- rep(social_dim_names[i],nrow(tbl))
  if(i==1){
    df_average_social <- tbl
  } else {
    df_average_social <- rbind(df_average_social,tbl)
  }
}

# Calculate range of ICCs over all permutations
df_single_emotion$icc_diff <- df_single_emotion$dim_icc_ub - df_single_emotion$dim_icc_lb
df_average_emotion$icc_diff <- df_average_emotion$dim_icc_ub - df_average_emotion$dim_icc_lb
df_single_social$icc_diff <- df_single_social$dim_icc_ub - df_single_social$dim_icc_lb
df_average_social$icc_diff <- df_average_social$dim_icc_ub - df_average_social$dim_icc_lb

# Number of subjects
df_single_emotion$subjects <- df_single_emotion$X+1
df_average_emotion$subjects <- df_average_emotion$X+1
df_single_social$subjects <- df_single_social$X+1
df_average_social$subjects <- df_average_social$X+1

# Plot only up to 20 subjects for emotion data
df_single_emotion <- df_single_emotion[!(df_single_emotion$subjects>20),]
df_average_emotion <- df_average_emotion[!(df_average_emotion$subjects>20),]

# Plot only up to 9 subjects for social data
df_single_social <- df_single_social[!(df_single_social$subjects>9),]
df_average_social <- df_average_social[!(df_average_social$subjects>9),]

# Factor 
df_single_emotion$subjects <- as.factor(df_single_emotion$subjects)
df_average_emotion$subjects <- as.factor(df_average_emotion$subjects)
df_single_social$subjects <- as.factor(df_single_social$subjects)
df_average_social$subjects <- as.factor(df_average_social$subjects)

# Plot emotion single measures ICC stability  as violin plot
#pdf(file = "/path/visualizations/FigSI5_iccstability/materials/icc_stability_single_emotion.pdf",width = 6,height=5)  
ggplot(df_single_emotion,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

# Plot emotion average measures ICC stability  as violin plot
#pdf(file = "/path/visualizations/FigSI5_iccstability/materials/icc_stability_average_emotion.pdf",width = 6,height=5)  
ggplot(df_average_emotion,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

# Plot social single measures ICC stability  as violin plot
#pdf(file = "/path/visualizations/FigSI5_iccstability/icc_stability_single_social.pdf",width = 6,height=5)  
ggplot(df_single_social,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

# Plot social average measures ICC stability  as violin plot
#pdf(file = "/path/visualizations/FigSI5_iccstability/icc_stability_average_social.pdf",width = 6,height=5)  
ggplot(df_average_social,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()


##-------------------------------------------------------------------------------------------------------------------------------------------
## 6. 3D PCoA of the features (supplementary media)

library(rgl)
library(lessR)
library(RColorBrewer)
library(stringr)

# Load hierarchical clustering data
load("/path/clustering_analyses/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

# Distance measure for PCoA
distance <- 1-cor(data) #Pearson

# PCoA
fit <- cmdscale(distance,eig=TRUE,k=135)

# Cluster membership by cluster (Hard-coded cluster indices from consensus hierarchical clustering in Figure 4)
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

# Create colors
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Create 3 different 3D animations to cover all PCs from 1 to 9
n <- 0
for( i in seq(from=1,to=3)){
  
  # Dimensions to plot
  n<-n+1
  x <- fit$points[,n]
  n<-n+1
  y <- fit$points[,n]
  n<-n+1
  z <- fit$points[,n]
  
  # Create 3D animation
  rgl_init <- function(new.device = FALSE, bg = "white", width = 1080) { 
    if( new.device | rgl.cur() == 0 ) {
      rgl.open()
      par3d(windowRect = 50 + c( 0, 0, width, width ) )
      rgl.bg(color = bg)
    }
    rgl.clear(type = c("shapes", "bboxdeco"))
    rgl.viewpoint(theta = 15, phi = 20, zoom = 0.6)
  }
  
  rgl_init(new.device = T,bg = "white")
  rgl.spheres(x, y, z, r = 0.01, color = cc)
  texts3d(x+0.03, y+0.03, z+0.03, text = colnames(data), col= "black", cex = 0.7)
  
  # Add axes and scale appropriately
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), ylim, c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 0), zlim, color = "black")
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[1], 0, 0), c(0, ylim[1], 0), 
                c(0, 0, zlim[1]))
  rgl.points(axes, color = "black", size = 5)
  
  # Names of the axes
  if(i==1){
    txt <- c("Pleasant - Unpleasant","Submission - Dominance","Impulsive - Cognitive")
    pc0 <- 1
    pc1 <- 3
  }
  if(i==2){
    txt <- c("Careless - Conscientious","Introversion - Extraversion","Playfulness - Sexuality")
    pc0 <- 4
    pc1 <- 6
  }
  
  if(i==3){
    txt <- c("Alone - With others","Femininty - Masculinity","PC9")
    pc0 <- 7
    pc1 <- 9
  }
  rgl.texts(axes, text = txt, color = "black",
            adj = c(0.5, -0.8), size = 5)
  
  aspect3d(1,1,1)
  
  widget <- rglwidget(width = 1920,height = 1080)
  
  #Save html
  filename = paste("/path/visualizations/Supplementary_media/",pc0,"_pcoa",pc1,".html",sep="")
  htmltools::save_html(widget, filename)
}




