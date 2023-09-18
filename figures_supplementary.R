# Supplementary Figures / materials
#
#       1. Estimate sufficient number of raters by ICC stability (FigSI1)
#         a. Based on appraisals data 
#         b. Based on social data (within sample stability)
#       2. Bimodality analysis and plot (FigSI2)
#       3. Heatmap of social feature ratings (FigSI3)
#       4. Polar plots of main HC dimensions (FigSI4)
#       5. Representative clustering results of the two datasets (FigSI5)
#       6. 3D htmls of the PCoA loadings (Supplementary material)

##-----------------------------------------------------------------------------------------------------------------------------------------------------
## 1.a. ICC stability plot (Based on appraisals data) (Figure SI1)

library(ggplot2)
library(stringr)

# Dimensions
dims <- dir("/path/icc_analysis_appraisals/results") # Folder of appraisals data (not public)
dim_names <- str_replace_all(dims,"_"," ")

# Read data for single measures
l <- c()
for(i in seq(from=1,to=length(dims))){
  tbl <- read.csv(paste("/path/icc_analysis_appraisals/results",dims[i],"/results_single_measures.csv",sep=""))
  l <- rbind(l,nrow(tbl))
  dim <- rep(dim_names[i],nrow(tbl))
  tbl$dim <- dim
  if(i==1){
    df <- tbl
  } else {
    df <- rbind(df,tbl)
  }
}

df$ci_diff <- df$dim_ci_ub - df$dim_ci_lb
df$icc_diff <- df$dim_icc_ub - df$dim_icc_lb
df$subjects <- df$X+1

#Exclude >20 subs
idx <- df$subjects>20
df2 <- df[!idx,]
df2$subjects <- as.factor(df2$subjects)

# Plot ICC stability as violin plot
#pdf(file = "icc_analysis_appraisals/icc_stability_single.pdf",width = 6,height=5)  
ggplot(df2,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

# Read data for average measures
for(i in seq(from=1,to=length(dims))){
  tbl <- read.csv(paste("/path/icc_analysis_appraisals/results",dims[i],"/results_average_measures.csv",sep=""))
  dim <- rep(dim_names[i],nrow(tbl))
  tbl$dim <- dim
  if(i==1){
    df <- tbl
  } else {
    df <- rbind(df,tbl)
  }
}

df$ci_diff <- df$dim_ci_ub - df$dim_ci_lb
df$icc_diff <- df$dim_icc_ub - df$dim_icc_lb
df$subjects <- df$X+1

#Exclude >20 subs
idx <- df$subjects>20
df2 <- df[!idx,]
df2$subjects <- as.factor(df2$subjects)

# Plot ICC stability as violin plot
#pdf(file = "/path/icc_analysis_appraisals/icc_stability_average.pdf",width = 6,height=5)  
ggplot(df2,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (average measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
#dev.off()

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1.ab ICC stability plot (Based on social data) (Figure SI1)

library(ggplot2)
library(stringr)

# Dimensions
dims <- dir("path/icc_analysis/results")
dim_names <- str_replace_all(dims,"_"," ")

# Read data for single measures
for(i in seq(from=1,to=length(dims))){
  tbl <- read.csv(paste("/path/icc_analysis/results/",dims[i],"/results_single_measures.csv",sep=""))
  dim <- rep(dim_names[i],nrow(tbl))
  tbl$dim <- dim
  if(i==1){
    df <- tbl
  } else {
    df <- rbind(df,tbl)
  }
}

df$ci_diff <- df$dim_ci_ub - df$dim_ci_lb
df$icc_diff <- df$dim_icc_ub - df$dim_icc_lb
df$subjects <- df$X+1

#Exclude >9 subs
idx <- df$subjects>9
df2 <- df[!idx,]
df2$subjects <- as.factor(df2$subjects)

# Plot ICC stability as violin plot
pdf(file = "/path/icc_analysis/icc_stability_single.pdf",width = 6,height=5)  
ggplot(df2,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (single measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
dev.off()

# Read data for average measures
for(i in seq(from=1,to=length(dims))){
  tbl <- read.csv(paste("/path/icc_analysis/results/",dims[i],"/results_average_measures.csv",sep=""))
  dim <- rep(dim_names[i],nrow(tbl))
  tbl$dim <- dim
  if(i==1){
    df <- tbl
  } else {
    df <- rbind(df,tbl)
  }
}

df$ci_diff <- df$dim_ci_ub - df$dim_ci_lb
df$icc_diff <- df$dim_icc_ub - df$dim_icc_lb
df$subjects <- df$X+1

#Exclude >9 subs
idx <- df$subjects>9
df2 <- df[!idx,]
df2$subjects <- as.factor(df2$subjects)

# Plot ICC stability as violin plot
pdf(file = "/path/icc_analysis/icc_stability_average.pdf",width = 6,height=5)  
ggplot(df2,aes(x=subjects,y=icc_diff)) +
  geom_violin(fill = "black") +
  ylim(0,0.75) +
  xlab("Number of subjects") + 
  ylab("ICC difference (average measures)") + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = "NA", size=0.5))
dev.off()

##-----------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Calculate and plot bimodality for the average ratings (FigSI2)

# Catenate data from different videosets
data <- read.csv("/path/data/average_ratings/average_ratings_videoset_1.csv")
data <- rbind(data,read.csv("/path/data/average_ratings/average_ratings_videoset_2.csv"))
data <- rbind(data,read.csv("/path/data/average_ratings/average_ratings_videoset_3.csv"))
data <- rbind(data,read.csv("/path/data/average_ratings/average_ratings_videoset_4.csv"))
data <- rbind(data,read.csv("/path/data/average_ratings/average_ratings_videoset_5.csv"))
data <- rbind(data,read.csv("/path/data/average_ratings/average_ratings_videoset_6.csv"))

# Binomial coefficient and Hartigan's dip test
res2 <- mt_check_bimodality(data,methods=c("BC", "HDS_sim"))

hds <- res2$HDS_p_value
hds_sim <- t(res2$HDS_simulated_p_value)
bc <- t(res2$BC)

# Fit Gaussian unimodal or bimodal distributions and test with AIC
aic_diff <- c()
for(i in seq(from=1,to=ncol(data))){
  fit_unimodal <- Mclust(data[,i],G=1,model="V")
  fit_bimodal <- Mclust(data[,i],G=2,model="V")
  aic_diff[i] <- AIC(fit_unimodal) - AIC(fit_bimodal)
}

stat <- bc*(1-hds_sim)
stat <- as.data.frame(stat)
stat$bc <- bc[,1]
stat$hds_sim <- hds_sim[,1]
stat$aic <- tanh(scale(aic_diff))
names <- rownames(stat)
stat$name <- str_replace_all(names,"_"," ")
colnames(stat) <- c("composite","bc","hds_sim","aic","names")
rm("names")

sign <- rep(" ",136)
sign[stat$bc>0.555 & stat$hds_sim<0.05] <- "**"
sign[stat$bc>0.555 & !stat$hds_sim<0.05] <- "*"

#png(file = "/path/visualizations/FigSIXX_bimodality/bmodality_bars_bc.png",height = 7500,width = 5000)
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
#dev.off()

##-----------------------------------------------------------------------------------------------------------------------------------------------------
## 3. Superheat population average ratings to visualize the variation in the ratings (FigSI3)

library(superheat)

# Load and order data
load("/path/analyses/analysis1a.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

#pdf("/path/visualizations/FigSIXX_superheat_ratings/superheat.pdf",height = 14,width = 16)
par(mar=c(2,10,2,2))
superheat(t(data),
          grid.hline = FALSE,
          force.left.label = TRUE,
          left.label.size = 0.3,
          left.label.text.size = 2.5,
          left.label.text.alignment = "right",
          #heat.col.scheme = "red",
          left.label.col = "white",
          #bottom.label.col = "white",
          heat.pal = c("white","#F4A582","#B2182B"),
          #heat.pal = "red",
          #heat.pal.values = c(0,0.5,1),
          heat.lim = c(0,100),
          legend = T)
#dev.off()

##---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. Polarplot main HC dimension based on their PCoA loadings (FigSI4)

library(ggplot2)
library(lessR)
library(corrplot)
library(tsne)

# Load hierarchical clustering data
load("/path/analyses/analysis1a.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

dims <- colnames(data)
dims <- str_replace_all(dims,"_"," ")

# Distance measure for PCoA
distance <- 1-cor(data) #Pearson

# PCoA
fit <- cmdscale(distance,eig=TRUE,k=135) # Only 56 first components have positive eigenvalues
loadings <- fit$points

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


names <- colnames(data)
df <- data.frame(names,loadings)
df$group <- factor(cc,levels = unique(cc))

# No clusters
idx <- (df$group==1 | df$group==2 | df$group==4 | df$group==5 | df$group==7 | df$group==12 | df$group==15 | df$group==17 | df$group==24 | df$group==29)
df1 <- df[!idx,]
names1 <- names[!idx]
idx <- unique(df1$group)

# Load cluster significance tables 
p_cluster <- read.csv("/path/pcoa/pcoa_cluster_loadings_pval.csv")
sign001_cluster <- read.csv("/path/pcoa/pcoa_cluster_loadings_pval_sing001.csv")
signbonf_cluster <- read.csv("/path/pcoa/pcoa_cluster_loadings_pval_signbonf.csv")
colnames(p_cluster) <-c("PC","Introversion","Feeding","Excecutive control","Comptence & ambition","Affection","Prosociality & pleasant feelings","Extraversion & humor","Masculinity","Verbal communication","Sympathetic activity","Gesturing","Nonverbal communication","Goal oriented behaviour","Antisociality","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical interaction")

#Create asterisk table
# Create asterisk matrix
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


pcoa <- as.factor(c("Pleasant - Unpleasant","Submission - Dominance","Impulsive - Cognitive","Careless - Conscientious","Introversion - Extraversion","Friendship - Sexuality","Alone - With others","Femininity - Masculinity"))
pcoa_neg <- c("Pleasant","Submission","Cognitive","Conscientious","Extraversion","Friendship","Alone","Masculinity")
pcoa_pos <- c("Unpleasant","Dominance","Impulsive","Careless","Introversion","Sexuality","With others","Femininity")
clusters <- c("Introversion","Feeding","Excecutive control","Comptence & ambition","Affection","Prosociality & pleasant feelings","Extraversion & humor","Masculinity","Verbal communication","Sympathetic activity","Gesturing","Nonverbal communication","Goal oriented behaviour","Antisociality","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical interaction")
d <- as.data.frame(pcoa)

cluster_loadings <- c()
require(lattice)
for(i in seq(1:20)){
  d$loadings <- colMeans(as.matrix(df1[df1$group==idx[i],2:9]))
  cluster_loadings <- rbind(cluster_loadings,d$loadings)
  d$signi <- sign[,i] 
  
  # Polar plot (direction colorcoded)
  #pdf(file = paste("/path/visualizations/Fig4_hc_pcoa_comparison/polars_noname/",clusters[i],".pdf",sep = ""),width = 10,height=10)
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
# 5. Representative clustering of the two datasets (movei clips / full movie) (FigSI5)

library(corrplot)

load("/path/analyses/analysis4.RData")
m_movie <- m_trim

load("/path/analyses/analysis5.RData")
m_megaperception <- m_trim

# Choose k for megaperception
#pdf("/path/visualizations/Fig5_clustervalidation/materials/consensus_megaperception_k16.pdf",width = 22,height = 20)
corrplot(m_megaperception, addrect = 16, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),)
#dev.off()

# Choose k for movie
#pdf("/path/visualizations/Fig5_clustervalidation/materials/consensus_movie_k19.pdf",width = 22,height = 20)
corrplot(m_movie, addrect = 19, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
#dev.off()

##-------------------------------------------------------------------------------------------------------------------------------------------
## 7. 3D PCoA of the features (supplementary material)

library(rgl)
library(lessR)
library(RColorBrewer)
library(stringr)

# Load hierarchical clustering data
load("/path/analyses/analysis1a.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]


distance <- 1-cor(data) #Pearson

fit <- cmdscale(distance,eig=TRUE,k=135)
x <- fit$points[,4]
y <- fit$points[,5]
z <- fit$points[,6]

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

# Create colors
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

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

# Add axes
lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}

xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
rgl.lines(xlim, c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), ylim, c(0, 0), color = "black")
rgl.lines(c(0, 0), c(0, 0), zlim, color = "black")

# Add a point at the end of each axes to specify the direction
axes <- rbind(c(xlim[1], 0, 0), c(0, ylim[1], 0), 
              c(0, 0, zlim[1]))
rgl.points(axes, color = "black", size = 5)

rgl.texts(axes, text = c("Careless - Conscientious","Introversion - Extraversion","Playfulness - Sexuality"), color = "black",
          adj = c(0.5, -0.8), size = 5)

aspect3d(1,1,1)

widget <- rglwidget(width = 1920,height = 1080)

#Save html
filename = "/path/pcoa/pcoa_3d_html/pcoa4_pcoa6.html"
htmltools::save_html(widget, filename)


