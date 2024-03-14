# Calculate single measures ICC, Cronbach's alpha, SEM and feature occurrence rate for movie frame experiment

# Severi Santavirta 26.2.2024

library(psych)
library(ltm)

# Catenate data from different videosets
data <- read.csv("/scratch/severi/projects/megaperception/data/data_git/data_avg_movieframes/average_ratings_videoset_1.csv")
data <- rbind(data,read.csv("/path/data_avg_movieframes/average_ratings_videoset_2.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieframes/average_ratings_videoset_3.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieframes/average_ratings_videoset_4.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieframes/average_ratings_videoset_5.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieframes/average_ratings_videoset_6.csv"))

# Clean feature names
dims <- colnames(data)
dims <- str_replace_all(dims,"_"," ")
dims <- str_replace_all(dims,"[.]","-")
dim_read <- colnames(data)
dim_read[92] <- "Exerting_self-control"

# Make a dataframe for both ICC methods and save for later use
icc2 <- c()
icc2_lb <- c()
icc2_ub <- c()
values <- c()
values_over <- c()
occurance <- c()
names <- c()
mu <- c()
mu_over <- c()
se <- c()
alpha <- c()
datapath <- "/path/data_csv_movieframes"
for(i in seq(1:length(data))){
  values <- c(values,data[,i])
  names <- c(names,rep(dims[i],468))
  mu <- c(mu,rep(mean(data[,i],na.rm=TRUE),468))
  
  v1 <- read.csv(paste(datapath,"/",dim_read[i],"_1.csv",sep = ""))
  v2 <- read.csv(paste(datapath,"/",dim_read[i],"_2.csv",sep = ""))
  v3 <- read.csv(paste(datapath,"/",dim_read[i],"_3.csv",sep = ""))
  v4 <- read.csv(paste(datapath,"/",dim_read[i],"_4.csv",sep = ""))
  v5 <- read.csv(paste(datapath,"/",dim_read[i],"_5.csv",sep = ""))
  v6 <- read.csv(paste(datapath,"/",dim_read[i],"_6.csv",sep = ""))
  
  # Calculate Cronbach's alpha
  alphaV1 <- cronbach.alpha(v1, standardized = FALSE, CI = FALSE)
  alphaV2 <- cronbach.alpha(v2, standardized = FALSE, CI = FALSE)
  alphaV3 <- cronbach.alpha(v3, standardized = FALSE, CI = FALSE)
  alphaV4 <- cronbach.alpha(v4, standardized = FALSE, CI = FALSE)
  alphaV5 <- cronbach.alpha(v5, standardized = FALSE, CI = FALSE)
  alphaV6 <- cronbach.alpha(v6, standardized = FALSE, CI = FALSE)
  alpha <- c(alpha,rep(mean(c(alphaV1$alpha,alphaV2$alpha,alphaV3$alpha,alphaV4$alpha,alphaV5$alpha,alphaV6$alpha)),468))
  
  # Occurrence rate
  occ <- occ_rate(v1)+occ_rate(v2)+occ_rate(v3)+occ_rate(v4)+occ_rate(v5)+occ_rate(v6)
  occurance <- c(occurance,rep(occ,468))
  
  # Average values
  over_rate <- c(over(v1),over(v2),over(v3),over(v4),over(v5),over(v6))
  values_over <- c(values_over,over_rate)
  over_rate <- over_rate[over_rate>0]
  mu_over <- c(mu_over,rep(mean(over_rate,na.rm=TRUE),468))
  
  # SEM
  se <- c(se,rep(mean(c(standardError(v1),standardError(v2),standardError(v3),standardError(v4),standardError(v5),standardError(v6))),468))
  
  # Create a combined design matrix for ICC (absolute) calculation over all videosets
  m <- matrix(NA,468,sum(length(v1),length(v2),length(v3),length(v4),length(v5),length(v6)))
  m[1:78,1:length(v1)] <- as.matrix(v1)
  m[79:156,(length(v1)+1):sum(length(v1),length(v2))] <- as.matrix(v2)
  m[157:234,(sum(length(v1),length(v2))+1):sum(length(v1),length(v2),length(v3))] <- as.matrix(v3)
  m[235:312,(sum(length(v1),length(v2),length(v3))+1):sum(length(v1),length(v2),length(v3),length(v4))] <- as.matrix(v4)
  m[313:390,(sum(length(v1),length(v2),length(v3),length(v4))+1):sum(length(v1),length(v2),length(v3),length(v4),length(v5))] <- as.matrix(v5)
  m[391:468,(sum(length(v1),length(v2),length(v3),length(v4),length(v5))+1):sum(length(v1),length(v2),length(v3),length(v4),length(v5),length(v6))] <- as.matrix(v6)
  
  # Calculate ICCs
  res <- ICC(m,alpha=.05)
  
  # Single measures ICCs
  icc2 <- c(icc2,rep(res$results$ICC[2],468))
  icc2_lb <-c(icc2_lb,rep(res$results$`lower bound`[2],468))
  icc2_ub <-c(icc2_ub,rep(res$results$`upper bound`[2],468))
  
}

# Save a ICC data frame with raw rating values for plotting
df_icc2 <- data.frame(values,values_over,occurance,mu,se,mu_over,icc2,icc2_lb,icc2_ub,alpha)
df_icc2$occurance_d <- occurance/468
df_icc2$names <- factor(names,levels = unique(names))
#write.csv(df_icc2,file = "/path/consistency_measures_for_plotting_movieframe.csv")

# Save also summary values for each video clip
ind <- seq(from=1,to=63648,by=468)
tbl_icc2 <- df_icc2[ind,seq(from=3,to=12)]
colnames(tbl_icc2) <- c("Occurrence","Mean values","SE","Mean values over 5","ICC","ICC_lb","ICC_ub","Alpha","Occurrence_d","Dimension")
#write.csv(tbl_icc2,file = "/path/consistency_measures_movieframe.csv")


