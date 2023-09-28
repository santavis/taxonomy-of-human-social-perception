# Figures

#   1. Ridgeplot average feature rating distributions (Figure 1)
#   2. ICC and occurrence rate figure (Figure 2)
#   3. Barplot PCoA loadings for main dimensions (Figure 3)
#   4. Correlation matrix of clusters (Figure 4)
#   5. Tsne plot (Figure 5)
#   6. Cluster generalization figure (Figure 6)
#        a) PCoA generalizatin plot (6a)
#        b) Clustering generalization (6b & 6c)

##---------------------------------------------------------------------------------------------------------------------
# 1. Ridgeplot average feature rarting distributions (Figure 1)

library(dplyr)
library(stringr)
library(ggridges)
library(ggplot2)
library(viridis)
library(psych)

# Load rating data 
df <- read.csv("/path/single_measures_icc_raw_values.csv")

# Arrange by mean rating
df <- arrange(df,mu_over)
df <- mutate(df,names=factor(names, levels=unique(names)))

# Delete rows where values_over = 0 (plot distributions only from non-zero values)
dff <- df[df$values_over>0,]

# Plot in three folds
l <- split(dff,f=dff$names)
require(lattice)
for(i in seq(from=1,to=3)) {
  if(i==1){
    df <- do.call(rbind.data.frame, l[1:45])
  }
  if(i==2){
    df <- do.call(rbind.data.frame, l[46:90])
  }
  if(i==3){
    df <- do.call(rbind.data.frame, l[91:136])
  }
  
  # Ridgeline plot with feature labels
  #pdf(file = paste("/path/visualizations/Fig1_ridgeline/fold",i,"_ridgeplot_names.pdf",sep=""),width = 4,height=12)  
  p <- ggplot(df, aes(x = values_over, y = names, fill = stat(x))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, bandwidth = 4,show.legend = F) +
    scale_fill_viridis_c(name = "Rating",option = "B") +
    theme(panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(size = 12))
  print(p)
  #dev.off()
  
  # Ridgeline plot without feature labels
  #pdf(file = paste("/path/visualizations/Fig1_ridgeline/fold",i,"_ridgeplot.pdf",sep=""),width = 4,height=12)
  p <- ggplot(df, aes(x = values_over, y = names, fill = stat(x))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, bandwidth = 4,show.legend = F) +
    scale_fill_viridis_c(name = "Rating",option = "B") +
    theme(panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank())
  print(p)
  #dev.off()
}

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## 2. Plot sorted occurance rate and ICC (figure 2)

library(dplyr)
library(stringr)
library(ggridges)
library(ggplot2)
library(viridis)
library(psych)

# Load ICC data 
df <- read.csv("/path/single_measures_icc_raw_values.csv")

#Sort by ICC and barplot
idx <- seq(from=1, to=31591, by=234)
df1 <- df[idx,]
df1 <- df1[order(df1$icc2),]
df1$names <- factor(df1$names,levels = unique(df1$names))

# Color ICC
cc <- rep(1,136)
cc[(df1$icc2>0.2 & df1$icc2<0.4)] <- 2
cc[(df1$icc2>0.4 & df1$icc2<0.6)] <- 3
cc[(df1$icc2>0.6 & df1$icc2<1)] <- 4
df1$group <- factor(cc,levels = unique(cc))

# Barplot ICC
#png("/path/visualizations/Fig2_iccbars/icc.png",width = 10000,height = 3000)
ggplot(df1,aes(x=names,y=icc2,width=1,fill = group)) +
  geom_bar(stat="identity") +
  #geom_linerange(aes(x=names,y=icc2, ymin=icc2_lb, ymax=icc2_ub)) +
  ylim(0,0.8) +
  scale_fill_manual(values = c("#D6604D","#F4A582","#92C5DE","#4393C3")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
#dev.off()

## Plot occurrence rate similarly
#pdf("/path/visualizations/Fig2_iccbars/occ.pdf",width = 10,height = 1)
ggplot(df1,aes(x=names,y=occurance_d)) +
  geom_point(size = 1) +
  #ylim(0,1) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),limits = c(0,1)) +
  #scale_fill_manual(values = c("#D6604D","#F4A582","#92C5DE","#4393C3")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
#dev.off()

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## 3. Barplot PCoA loadings with differend N of features for each components (Figure 3)

library(ggplot2)
library(stringr)
library(lessR)

# Load hierarchical clustering data
load("/path/clustering_analyses/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

dims <- colnames(data)
dims <- str_replace_all(dims,"_"," ")

# Distance measure for PCoA
distance <- 1-cor(data) #Pearson

# PCoA
fit <- cmdscale(distance,eig=TRUE,k=135) # Only 56 first components have positive eigenvalues
loadings <- fit$points

# Compute variance explained
eig <- fit$eig
exp <- eig/sum(eig[eig>0]) # Scale eigenvalues to the sum of positive eigenvalues -> variance explained
exp_cum <- cumsum(exp)

# Create a data frame
df <- data.frame(dims,loadings[,1:8])
df <- df[order(df$dims),]

# Load significance tables 
sign001 <- read.csv("/path/pcoa/pcoa_loadings_pval_sign001.csv")
signbonf <- read.csv("/path/pcoa/pcoa_loadings_pval_signbonf.csv")
sign001 <- sign001[order(sign001$X),]
signbonf <- signbonf[order(signbonf$X),]

# Give significance color for each loading
sign <- matrix("",nrow=136,ncol=8,byrow=TRUE)
sign[as.matrix(df[,2:9])>0] <- "#FDDBC7"
sign[as.matrix(df[,2:9])<0] <- "#D1E5F0"
sign[as.matrix(sign001[,2:9]) & as.matrix(df[,2:9])>0] <- "#D6604D"
sign[as.matrix(sign001[,2:9]) & as.matrix(df[,2:9])<0] <- "#4393C3"
sign[as.matrix(signbonf[,2:9]) & as.matrix(df[,2:9])>0] <- "#B2182B"
sign[as.matrix(signbonf[,2:9]) & as.matrix(df[,2:9])<0] <- "#2166AC"
df <- cbind(df,as.data.frame(sign))
rm("sign")

# How many bars to plot per dimension
n <- c(20,10,10,10,10,10,4,3)
h <- c(3000,1500,1500,1500,1500,1500,600,450)
require(lattice)
for(i in seq(from=1,to=8)){ # Plot only significant components (first 8)
  
  #Create a new daata frame from the ith components data
  df1 <- df[,c(1,i+1,i+1+8)]
  df1 <- df1[order(df1[,2]),] # Sort by loading
  idx <- c(seq(from=1,to=n[i]),seq(from=(136-n[i]+1),to=136)) # Select which bars to plot
  df1 <- df1[idx,] # Select which bars to plot
  df1$dims <- factor(df1$dims,levels = df1$dims)
  colnames(df1) <- c("d","load","sign")
  
  # Plot component names
  png(file = paste("/path/visualizations/Fig4_pcoa/barplots/",i,"_component_name.png",sep = ""),width = 1000,height=h[i])
  p <- ggplot(df1,aes(x=d,y=load,fill=sign)) +
    geom_bar(stat = "identity") +
    scale_fill_identity() +
    ylim(-1,1) +
    coord_flip() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 16))
  print(p)
  dev.off()
  
  # Plot component loadings
  png(file = paste("/path/visualizations/Fig4_pcoa/barplots/",i,"_component.png",sep = ""),width = 1000,height=h[i])
  p <- ggplot(df1,aes(x=d,y=load,fill=sign)) +
    geom_bar(stat = "identity") +
    scale_fill_identity() +
    ylim(-1,1) +
    coord_flip() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 16))
  print(p)
  dev.off()
}


##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 4. Results from the hierarchical clustering analysis (analysis1) (Figure 4)

library(corrplot)
library(lessR)

load("/path/clustering_analyses/analysis1.RData")
ordered <- corReorder(m_trim, order = "hclust", hclust_type = "average")
data <- data[,colnames(ordered)]

# Plot analysis 1 correlation matrix ordered to match with the consensus matrix
cormat <- cor(data)

#cormat
#pdf(file = "/path/visualizations/Fig4_clustermatrix/materials/cormat_clips.pdf",width = 22,height=20)
corrplot(cormat,tl.col = "black",method = "square", col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
#dev.off()

#cormat lower
#pdf(file = "/path/visualizations/Fig4_clustermatrix/materials/cormat_lower_clips.pdf",width = 22,height=20)
corrplot(cormat,tl.col = "black",method = "square", col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),type = 'lower',tl.pos = "n")
#dev.off()

#consensus upper 
#pdf(file = "/path/visualizations/Fig4_clustermatrix/materials/consensus_upper_clips.pdf",width = 22,height=20)
corrplot(ordered,tl.col = "black",method = "square",col.lim = c(0, 1), col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),type = 'upper',tl.pos = "n")
#dev.off()

#Consensus with rectangles
#pdf(file = "/path/visualizations/Fig4_clustermatrix/materials/consensus_30rect_clips.pdf",width = 22,height=20)
corrplot(ordered,order = "hclust", hclust.method = "average", tl.col = "black",method = "square",col.lim = c(0, 1), col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),tl.pos = "n",addrect = 30)
#dev.off()


##---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5. Tsne plot (Figure 5).

library(ggplot2)
library(lessR)
library(corrplot)
library(tsne)
library(pals)

# Load hierarchical clustering data
load("/path/clustering_analyses/analysis1.RData")
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

# Make data frame
names <- colnames(data)
df <- data.frame(names,loadings[,1:8])
df$group <- factor(cc,levels = unique(cc))

# Exclude unclustered features
idx <- (df$group==1 | df$group==2 | df$group==4 | df$group==5 | df$group==7 | df$group==12 | df$group==15 | df$group==17 | df$group==24 | df$group==29)
df1 <- df[!idx,]
names1 <- names[!idx]
idx <- unique(df1$group)

pcoa <- as.factor(c("Pleasant - Unpleasant","Submission - Dominance","Impulsive - Cognitive","Careless - Conscientious","Introversion - Extraversion","Friendship - Sexuality","Alone - With others","Femininity - Masculinity"))
pcoa_neg <- c("Pleasant","Submission","Cognitive","Conscientious","Extraversion","Friendship","Alone","Masculinity")
pcoa_pos <- c("Unpleasant","Dominance","Impulsive","Careless","Introversion","Sexuality","With others","Femininity")
clusters <- c("Introversion","Feeding","Excecutive control","Comptence & ambition","Affection","Prosociality & pleasant feelings","Extraversion & humor","Masculinity","Verbal communication","Sympathetic activity","Gesturing","Nonverbal communication","Goal oriented behaviour","Antisociality","Unpleasant feelings","Body movement","Submission","Sexuality","Femininity","Physical interaction")
d <- as.data.frame(pcoa)

# Plot cluster members in 2D based on TNSE of the loadings of first eigth PCs
tsn <- tsne(as.matrix(df1[2:9]), initial_dims = 8)
tsn <- data.frame(tsn)
tsn$gr <- factor(df1$group,levels = unique(df1$group))

#Try plotting cluster number over the point to make it easier to distinguish between the colours
n <- rep(0,126)
n[tsn$gr==3] <- 1
n[tsn$gr==6] <- 2
n[tsn$gr==8] <- 3
n[tsn$gr==9] <- 4
n[tsn$gr==10] <- 5
n[tsn$gr==11] <- 6
n[tsn$gr==13] <- 7
n[tsn$gr==14] <- 8
n[tsn$gr==16] <- 9
n[tsn$gr==18] <- 10
n[tsn$gr==19] <- 11
n[tsn$gr==20] <- 12
n[tsn$gr==21] <- 13
n[tsn$gr==22] <- 14
n[tsn$gr==23] <- 15
n[tsn$gr==25] <- 16
n[tsn$gr==26] <- 17
n[tsn$gr==27] <- 18
n[tsn$gr==28] <- 19
n[tsn$gr==30] <- 20
tsn$grnumber <- as.factor(n)

# Get colours to the points
pal <- alphabet(n = 20)
cv <- rep("",126)
cv[tsn$gr==3] <- pal[1]
cv[tsn$gr==6] <- pal[2]
cv[tsn$gr==8] <- pal[3]
cv[tsn$gr==9] <- pal[4]
cv[tsn$gr==10] <- pal[5]
cv[tsn$gr==11] <- pal[6]
cv[tsn$gr==13] <- pal[7]
cv[tsn$gr==14] <- pal[8]
cv[tsn$gr==16] <- pal[9]
cv[tsn$gr==18] <- pal[10]
cv[tsn$gr==19] <- pal[11]
cv[tsn$gr==20] <- pal[12]
cv[tsn$gr==21] <- pal[13]
cv[tsn$gr==22] <- pal[14]
cv[tsn$gr==23] <- pal[15]
cv[tsn$gr==25] <- pal[16]
cv[tsn$gr==26] <- pal[17]
cv[tsn$gr==27] <- pal[18]
cv[tsn$gr==28] <- pal[19]
cv[tsn$gr==30] <- pal[20]
tsn$cc <- cv

#pdf("/path/visualizations/Fig5_tsne/tsne_numbered_dots_new.pdf",width = 8,height = 8)
ggplot() +
  geom_point(aes(x=tsn$X1,y=tsn$X2,color=tsn$grnumber),size=7,show.legend = FALSE) +
  scale_color_manual(labels = clusters, values = unique(tsn$cc)) + 
  geom_text(aes(x=tsn$X1,y=tsn$X2, label = tsn$grnumber),size=4,fontface = "bold",color = "white") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text=element_text(size=15))
#dev.off()

##---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6a. PCoA generalization plot (Figure 7a) 
library(corrplot)
library(lessR)

# Load data
load("/path/clustering_analyses/analysis2.RData")
data_movie <- data
load("/path/clustering_analyses/analysis3.RData")
data_clips <- data
data_movie <- data_movie[,colnames(data_clips)]

# Distance measure
distance_clips <- 1-cor(data_clips) #Pearson
distance_movie <- 1-cor(data_movie) #Pearson

# PCoA with real data
fit_clips <- cmdscale(distance_clips,eig=TRUE,k=75)
fit_movie <- cmdscale(distance_movie,eig=TRUE,k=75)
loadings_clips <- fit_clips$points
loadings_movie <- fit_movie$points
weights_movie <- fit_movie$eig
weights_clips <- fit_clips$eig
var_exp_movie <- weights_movie[weights_movie>0]/sum(weights_movie[weights_movie>0])
var_exp_clips <- weights_clips[weights_clips>0]/sum(weights_clips[weights_clips>0])

# Take first 7 dimensions from clips and 10 dimensions from movie (signficant components)
loadings7_clips <- loadings_clips[,1:7]
cats <- c("Clips PC1","Clips PC2","Clips PC3","Clips PC4","Clips PC5","Clips PC6","Clips PC7") 
colnames(loadings7_clips) <- cats

loadings10_movie <- loadings_movie[,1:10]
loadings10_movie[,1] <- loadings10_movie[,1]*-1 # Flip weight signs to avoid high negative correlations in the correlations matrix. The direction of the axis is insignificant
loadings10_movie[,5] <- loadings10_movie[,5]*-1 # Flip weight signs to avoid high negative correlations in the correlations matrix. The direction of the axis is insignificant
loadings10_movie[,9] <- loadings10_movie[,9]*-1 # Flip weight signs to avoid high negative correlations in the correlations matrix. The direction of the axis is insignificant
cats <- c("Movie PC1","Movie PC2","Movie PC3","Movie PC4","Movie PC5","Movie PC6","Movie PC7","Movie PC8","Movie PC9","Movie PC10")
colnames(loadings10_movie) <- cats

loadings <- cbind(loadings7_clips,loadings10_movie)
cormat <- corReorder(cor(loadings),order = "hclust")
loadings <- loadings[,colnames(cormat)]
cortest <- cor.mtest(loadings, conf.level = 0.95)

#pdf("/path/visualizations/Fig6_generalization/materials/pcoa_generalization.pdf",width = 8,height = 8)
corrplot(cormat, order = "hclust",hclust.method = "average", type = "lower", p.mat = cortest$p, insig='blank', tl.col = ifelse(grepl("Clips",colnames(cormat)),"#2166AC","#B2182B"), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),tl.cex = 1.5)
#dev.off()

##---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6b. Clustering generalization plot, alluvial & correlation matrices (Figure 6b & 6c)

library(corrplot)
library(ggplot2)
library(ggalluvial)
library(stringr)
library(lessR)

load("/path/clustering_analyses/analysis2.RData")
m_movie <- m_trim

load("/path/clustering_analyses/analysis3.RData")
m_clips <- m_trim

# Check for best number of clusters based on corrplots
b <- corrplot(m_clips, addrect = 16, order = "hclust", hclust.method = "average", tl.col = "black",method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
order_clips <- b$corr
clips_clusters <- factor(c(1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,4,4,4,5,6,6,6,6,7,7,7,7,7,8,8,8,8,9,10,11,11,11,11,12,12,13,13,13,13,13,13,13,13,13,14,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16),level = c(1:16))

c <- corrplot(m_movie, addrect = 19, order = "hclust", hclust.method = "average", tl.col = "black",method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
order_movie <- c$corr

movie_clusters <- factor(c(1,1,2,2,3,3,3,3,4,4,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,8,9,9,9,9,9,9,9,9,10,11,12,12,12,12,13,13,13,14,14,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,18,19,19,19),level = c(1:19))
movie_clusters_order <- c(1:76)
movie_clusters_order[movie_clusters==16] = 1
movie_clusters_order[movie_clusters==14] = 2
movie_clusters_order[movie_clusters==19] = 3
movie_clusters_order[movie_clusters==17] = 4
movie_clusters_order[movie_clusters==12] = 5
movie_clusters_order[movie_clusters==1] = 6
movie_clusters_order[movie_clusters==8] = 7
movie_clusters_order[movie_clusters==9] = 8
movie_clusters_order[movie_clusters==4] = 9
movie_clusters_order[movie_clusters==18] = 10
movie_clusters_order[movie_clusters==11] = 11
movie_clusters_order[movie_clusters==15] = 12
movie_clusters_order[movie_clusters==13] = 13
movie_clusters_order[movie_clusters==7] = 14
movie_clusters_order[movie_clusters==3] = 15
movie_clusters_order[movie_clusters==10] = 16
movie_clusters_order[movie_clusters==2] = 17
movie_clusters_order[movie_clusters==5] = 18
movie_clusters_order[movie_clusters==6] = 19
movie_clusters_order <- factor(movie_clusters_order,level = c(1:19))

#Flip horizontally and vertically to better match clips in alluvial
order_movie_flip <- apply(order_movie, 2, rev)
order_movie_flip <- apply(order_movie_flip, 1, rev)
corrplot(order_movie_flip,tl.col = "black",method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
movie_clusters_flip <- factor(c(1,1,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,6,6,7,7,7,8,8,8,8,9,10,11,11,11,11,11,11,11,11,12,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,16,16,17,17,17,17,18,18,19,19),level = c(1:19))

tmp <- data.frame(movie_clusters)
rownames(tmp) <- colnames(order_movie)
movie_clusters <- tmp[colnames(order_clips),]

tmp <- data.frame(movie_clusters_order)
rownames(tmp) <- colnames(order_movie)
movie_clusters_order <- tmp[colnames(order_clips),]

tmp <- data.frame(movie_clusters_flip)
rownames(tmp) <- colnames(order_movie_flip)
movie_clusters_flip <- tmp[colnames(order_clips),]

#Plot alluvial
dims <- colnames(order_clips)
alluvial_data <- data.frame(dims,clips_clusters,movie_clusters,movie_clusters_order,movie_clusters_flip)
alluvial_data$id <- factor(c(1:76),level=c(1:76))

remove(movie_clusters,movie_clusters_order,clips_clusters,movie_clusters_flip)
group = factor(c(1:35),level = c(1:35))

#pdf("/path/visualizations/Fig6_generalization/materials/validation_alluvial.pdf",width = 20,height = 20)
ggplot(alluvial_data,aes(y = rep(1,76),axis1 = clips_clusters, axis2 = movie_clusters_flip)) +
  scale_x_discrete(limits = c("clips", "Movie")) +
  geom_alluvium(aes(fill = clips_clusters),width = 0.3,) +
  #stat_alluvium(geom = "text", aes(label = dims)) +
  theme_void() +
  theme(legend.position = "none")
#dev.off()

#pdf("/path/visualizations/Fig6_generalization/materials/validation_alluvial_rectangles.pdf",width = 20,height = 20)
ggplot(alluvial_data,aes(y = rep(1,76),axis1 = clips_clusters, axis2 = movie_clusters_flip)) +
  geom_stratum(width = 0.3) +
  theme_void() +
  theme(legend.position = "none")
#dev.off()

# Choose k for clips
#pdf("/path/visualizations/Fig6_generalization/materials/consensus_clips_k16.pdf",width = 22,height = 20)
corrplot(m_clips, addrect = 16, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),)
#dev.off()

# Choose k for movie
#pdf("/path/visualizations/Fig6_generalization/materials/consensus_movie_k19.pdf",width = 22,height = 20)
corrplot(m_movie, addrect = 19, order = "hclust", hclust.method = "average", tl.col = "black", col.lim = c(0, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20))
#dev.off()

#pdf("/path/visualizations/Fig6_generalization/materials/cormat_movie.pdf",width = 22,height = 20)
corrplot(order_movie, tl.col = "black", col.lim = c(-1, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),type = "upper")
#dev.off()

#pdf("/path/visualizations/Fig6_generalization/materials/cormat_clips.pdf",width = 22,height = 20)
corrplot(order_clips, tl.col = "black", col.lim = c(-1, 1), method = "square",col=colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))(20),type = "lower")
#dev.off()





