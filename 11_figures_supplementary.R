# Supplementary figures / media
#
#       1. Bimodality analysis and plot (FigSI1)
#       2. Heatmap of social feature ratings (FigSI2)
#       3. Concordance analysis, polar plots of main HC dimensions (FigSI3)
#       4. Representative clustering results of the movie clip and movie frame data in Study 3 (FigSI4)
#       5. Estimate sufficient number of raters by ICC stability (FigSI5)
#       6. Demographical analysis of the social ratings (FigSI6)
#       7. 3D htmls of the PCoA loadings (Supplementary media)
#
# Severi Santavirta 7.5.2024

##-----------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Calculate and plot bimodality for the average ratings (FigSI1)

library(ggplot2)
library(mousetrap)
library(stringr)

# Catenate data from different videosets
data <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_1.csv")
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_2.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_3.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_4.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_5.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_6.csv"))

# Binomial coefficient and Hartigan's dip test
res <- mt_check_bimodality(data,methods=c("BC", "HDS_sim"))

# Create a data frame for plotting
stat <- data.frame(t(res$BC)[,1])
stat$hds_sim <- t(res$HDS_simulated_p_value)[,1]
stat$name <- str_replace_all(colnames(data),"_"," ")
colnames(stat) <- c("bc","hds_sim","names")

# Create significance asterisks
sign <- rep(" ",136)
sign[stat$bc>0.555 & stat$hds_sim<0.05] <- "**"
sign[stat$bc>0.555 & !stat$hds_sim<0.05] <- "*"

# Barplot sorted Bimodality coeffs witth significance 
#png(file = "/path/FigSI1_bimodality/FigSI1_bimodality_bars_bc.png",height = 7500,width = 5000)
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
load("/path/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

#pdf("/path/FigSI2_superheat_ratings/superheat.pdf",height = 14,width = 16)
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
library(tsne)

# Load hierarchical clustering data and correct names
load("/path/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]
dims <- colnames(data)
dims <- str_replace_all(dims,"_"," ")

# Distance measure for PCoA
distance <- 1-cor(data) #Pearson

# PCoA
fit <- cmdscale(distance,eig=TRUE,k=135)
loadings <- fit$points

# Cluster membership by cluster (Hard-coded cluster indices from consensus hierarchical clustering
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
p_cluster <- read.csv("/path/pcoa_cluster_loadings_pval.csv")
sign001_cluster <- read.csv("/path/pcoa_cluster_loadings_pval_sing001.csv")
signbonf_cluster <- read.csv("/path/pcoa_cluster_loadings_pval_signbonf.csv")
colnames(p_cluster) <-c("PC","Introversion","Feeding","Self-control","Achievement","Emotional affection","Pleasant feelings & prosociality","Extraversion & playfulness","Masculinity","Social engagement","Physical discomfort","Gesturing","Emotional expression","Motivation","Antisocial behaviour","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical affection")

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
pcoa <- as.factor(c("Pleasant - Unpleasant","Empathy - Dominance","Physical - Cognitive","Inactive - Proactive","Introversion - Extraversion","Playful - Sexual","Alone - Interaction","Femininity - Masculinity"))
pcoa_neg <- c("Pleasant","Empathy","Cognitive","Proactive","Extraversion","Playful","Alone","Masculinity")
pcoa_pos <- c("Unpleasant","Dominance","Physical","Inactive","Introversion","Sexual","Interaction","Femininity")
clusters <- c("Introversion","Feeding","Self-control","Achievement","Emotional affection","Pleasant feelings & prosociality","Extraversion & playfulness","Masculinity","Social engagement","Physical discomfort","Gesturing","Emotional expression","Motivation","Antisocial behaviour","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical affection")
d <- as.data.frame(pcoa)

# Polarplot main dimension loadings for each HC cluster
cluster_loadings <- c()
require(lattice)
for(i in seq(1:20)){
  d$loadings <- colMeans(as.matrix(df1[df1$group==idx[i],2:9]))
  cluster_loadings <- rbind(cluster_loadings,d$loadings)
  d$signi <- sign[,i] 
  
  # Polar plot (direction colorcoded)
  #pdf(file = paste("/path/FigSI3_polarplots/polars_noname/",clusters[i],".pdf",sep = ""),width = 10,height=10)
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
# 4. Representative clustering of the two datasets (movie clips / full movie) that serves as the ground for the alluvial diagram in Figure 9 (FigSI4)

library(corrplot)

load("/path/analysis2.RData")
m_movie <- m_trim

load("/path/analysis3.RData")
m_megaperception <- m_trim

# Choose k for megaperception
#pdf("/path/FigSI4_generalization_consensusmatrices/materials/consensus_megaperception_k16.pdf",width = 22,height = 20)
corrplot(m_megaperception, addrect = 16, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),)
#dev.off()

# Choose k for movie
#pdf("/path/FigSI4_generalization_consensusmatrices/materials/consensus_movie_k19.pdf",width = 22,height = 20)
corrplot(m_movie, addrect = 19, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
#dev.off()

##-----------------------------------------------------------------------------------------------------------------------------------------------------
## 5. ICC stability plot (Based on appraisals data) (Figure SI5)

library(ggplot2)
library(stringr)

# Emotion dimensions (data not publicly available)
emotion_dims <- dir("/path/icc_analysis_appraisals/results") # Path to appraisals data
emotion_dim_names <- str_replace_all(emotion_dims,"_"," ")

# Social dimensions
social_dims <- dir("/path/icc_analysis/results")
social_dim_names <- str_replace_all(social_dims,"_"," ")


# Read ICC data for emotions (data not publicly available)
for(i in seq(from=1,to=length(emotion_dims))){
  
  # Emotion single measures ICC
  tbl <- read.csv(paste("/path/icc_analysis_appraisals/results/",emotion_dims[i],"/results_single_measures.csv",sep=""))
  tbl$dim <- rep(emotion_dim_names[i],nrow(tbl))
  if(i==1){
    df_single_emotion <- tbl
  } else {
    df_single_emotion <- rbind(df_single_emotion,tbl)
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
}

# Calculate range of ICCs over all permutations
df_single_emotion$icc_diff <- df_single_emotion$dim_icc_ub - df_single_emotion$dim_icc_lb
df_single_social$icc_diff <- df_single_social$dim_icc_ub - df_single_social$dim_icc_lb

# Number of subjects
df_single_emotion$subjects <- df_single_emotion$X+1
df_single_social$subjects <- df_single_social$X+1

# Plot only up to 20 subjects for emotion data
df_single_emotion <- df_single_emotion[!(df_single_emotion$subjects>20),]

# Plot only up to 9 subjects for social data
df_single_social <- df_single_social[!(df_single_social$subjects>9),]

# Factor 
df_single_emotion$subjects <- as.factor(df_single_emotion$subjects)
df_single_social$subjects <- as.factor(df_single_social$subjects)

# Plot emotion single measures ICC stability  as violin plot
#pdf(file = "/path/FigSI5_iccstability/materials/icc_stability_single_appraisals.pdf",width = 6,height=5)  
ggplot(df_single_emotion,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

# Plot social single measures ICC stability  as violin plot
#pdf(file = "/path/FigSI5_iccstability/icc_stability_single_social.pdf",width = 6,height=5)  
ggplot(df_single_social,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 7. Correlate each participants responses to the mean of all others (ethnicity % age related features).
# For sex we have balanced sample so we compare whether the males ratings correlate more with the means of other males than with the means of females
# The goal is to identify if there are major differences in the perception between different demographic groups (FigSI6)
# The demographical data is not shared publicly since it requires to link public Prolific ID with the ratings.

library(stringr)
library(ggplot2)

# Input data
demo <- read.csv("/path//participants_analyzed.csv") # Data is not shared
data_path <- "/path/data_csv"

# Read feature names from clustering data
load("/path/analysis1.RData")
dims <- colnames(data)
dims <- str_replace_all(dims," ","_")


# Loop thorough features
participant <- c()
feature <- c()
sex <- c()
age <- c()
race <- c()
corr <- c()
corr_same_sex <- c()
corr_other_sex <- c()
n <- 0
for(feat in seq(from=1,to=length(dims))){
  
  # Loop through all video sets
  for(set in seq(from=1,to=6)){
    
    # Read data for the specific feature and video set
    data <- read.csv(paste(data_path,"/",dims[feat],"_",set,".csv",sep = ""))
    
    # Read participant names
    subj <- colnames(data)
    subj <- substr(subj,2,nchar(subj))
    
    # Get the sex information of all subjects
    sex_set <- c()
    for(s in seq(from=1,to=length(subj))){
      idx <- which(demo$ParticipantId==subj[s])
      sex_set[s] <- demo$Sex[idx]
    }
    
    # Loop thorough subjects
    for(s in seq(from=1,to=length(subj))){
      
      # Store information
      n <- n+1
      participant[n] <- subj[s]
      feature[n] <- dims[feat]
      
      # Find demographics of this participant
      idx <- which(demo$ParticipantId==subj[s])
      sex[n] <- demo$Sex[idx]
      age[n] <- demo$Age[idx]
      race[n] <- demo$EthnicitySimplified[idx]
      
      # Calculate the correlation of the ratings between all others
      corr[n] <- cor(data[,s],rowMeans(data[,-s]))
      
      # Correlation between other sex
      if(sex[n]=="Male"){
        if(sum(sex_set=="Female")>1){
          corr_other_sex[n] <- cor(data[,s],rowMeans(data[,sex_set=="Female"]))
        }
        if(sum(sex_set=="Female")==1){
          corr_other_sex[n] <- cor(data[,s],data[,which(sex_set=="Female")])
        }
        if(sum(sex_set=="Female")==0){
          corr_other_sex[n] <- NA
        }
      }
      if(sex[n]=="Female"){
        if(sum(sex_set=="Male")>1){
          corr_other_sex[n] <- cor(data[,s],rowMeans(data[,sex_set=="Male"]))
        }
        if(sum(sex_set=="Male")==1){
          corr_other_sex[n] <- cor(data[,s],data[,which(sex_set=="Male")])
        }
        if(sum(sex_set=="Male")==0){
          corr_other_sex[n] <- NA
        }
      }
      if(sex[n]=="DATA_EXPIRED"){
        corr_other_sex[n] <- NA
      }
      
      # Correlation between same sex
      data_same_sex <- data[,-s]
      sex_set_s <- sex_set[-s]
      if(sex[n]=="Male"){
        if(sum(sex_set_s=="Male")>1){
          corr_same_sex[n] <- cor(data[,s],rowMeans(data_same_sex[,sex_set_s=="Male"]))
        }
        if(sum(sex_set_s=="Male")==1){
          corr_same_sex[n] <- cor(data[,s],data_same_sex[,which(sex_set_s=="Male")])
        }
        if(sum(sex_set_s=="Male")==0){
          corr_same_sex[n] <- NA
        }
      }
      if(sex[n]=="Female"){
        if(sum(sex_set_s=="Female")>1){
          corr_same_sex[n] <- cor(data[,s],rowMeans(data_same_sex[,sex_set_s=="Female"]))
        }
        if(sum(sex_set_s=="Female")==1){
          corr_same_sex[n] <- cor(data[,s],data_same_sex[,which(sex_set_s=="Female")])
        }
        if(sum(sex_set_s=="Female")==0){
          corr_same_sex[n] <- NA
        }
      }
      if(sex[n]=="DATA_EXPIRED"){
        corr_same_sex[n] <- NA
      }
    }
  }
}

output <- data.frame(corr)
colnames(output) <- c("r")
output$r_same_sex <- corr_same_sex
output$r_other_sex <- corr_other_sex
output$feature <- feature
output$sex <- sex
output$age <- age
output$race <- race

# Remove rows where we do not have complete data
idx_age <- which(output$age==0 | is.nan(output$age))
output <- output[-idx_age,]
idx_race <- which(output$race=="DATA_EXPIRED")
output <- output[-idx_race,]
idx_r <- which(is.na(output$r))
output <- output[-idx_r,]
idx_r_same_sex <- which(is.na(output$r_same_sex))
output <- output[-idx_r_same_sex,]

# Association between ratings similarity and age
age_effect <- cor(output$r,output$age)
#pdf('/path/age_scatter.pdf',height = 3,width = 3)
ggplot(output, aes(x = age, y = r)) +
  geom_point(size=0.5) +
  geom_smooth(method = "loess", se = FALSE) +  # Adding local regression line
  labs(y = "Rating similarity with others (r)", x = "Age") +
  theme_minimal()
#dev.off()

# Association between rating similarity and race
means_race <- aggregate(r ~ race, data = output, FUN = mean)
#pdf('/path/race_violin.pdf',height = 3,width = 3)
ggplot(output, aes(x = race, y = r)) +
  geom_violin(fill = "#F4A582",show.legend = FALSE) +
  geom_point(data = means_race, aes(x = race, y = r), color = "black", size = 3) +  # Overlay means as points
  labs(y = "Rating similarity with others (r)", x = "Ethnicity") +
  theme_minimal()
#dev.off()

output_sex <- c(output$r_same_sex,output$r_other_sex)
output_sex <- data.frame(output_sex)
colnames(output_sex) <- c("r")
output_sex$r_type <- c(rep("Same sex",nrow(output)),rep("Opposite sex",nrow(output)))

# Association between rating similarity and sex
means_sex <- aggregate(r ~ r_type, data = output_sex, FUN = mean)
#pdf('/path/sex_violin.pdf',height = 3,width = 3)
ggplot(output_sex, aes(x = r_type, y = r)) +
  geom_violin(fill = "#92C5DE",show.legend = FALSE) +
  geom_point(data = means_sex, aes(x = r_type, y = r), color = "black", size = 3) +  # Overlay means as points
  labs(y = "Rating similarity with other sex (r)", x = "Sex") +
  theme_minimal()
#dev.off()

##-------------------------------------------------------------------------------------------------------------------------------------------
## 7. 3D PCoA of the features (supplementary media)

library(rgl)
library(lessR)
library(RColorBrewer)
library(stringr)

# Load hierarchical clustering data
load("/path/analysis1.RData")
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
    txt <- c("Pleasant - Unpleasant","Empathy - Dominance","Physical - Cognitive")
    pc0 <- 1
    pc1 <- 3
  }
  if(i==2){
    txt <- c("Inactive - Proactive","Introversion - Extraversion","Playful - Sexual")
    pc0 <- 4
    pc1 <- 6
  }
  
  if(i==3){
    txt <- c("Alone - Interaction","Femininity - Masculinity","PC9")
    pc0 <- 7
    pc1 <- 9
  }
  rgl.texts(axes, text = txt, color = "black",
            adj = c(0.5, -0.8), size = 5)
  
  aspect3d(1,1,1)
  
  widget <- rglwidget(width = 1920,height = 1080)
  
  #Save html
  filename = paste("/path/Supplementary_media/pcoa",pc0,"_pcoa",pc1,".html",sep="")
  htmltools::save_html(widget, filename)
}



