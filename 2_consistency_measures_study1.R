# Calculate single measures ICC, Cronbach's alpha, SEM and feature occurrence rate for the primary movie clip data (Study 1)

# Severi Santavirta 26.2.2024

library(psych)
library(ltm)

# Calculate occurance rate
occ_rate <- function(v){
  mu <- rowMeans(v)
  sem <- apply(v,1,sd)/sqrt(ncol(v))
  low <- mu-sem
  occ <- sum(low>5)
  return(occ)
}

# Return mu_over
over <- function(v){
  mu <- rowMeans(v)
  sem <- apply(v,1,sd)/sqrt(ncol(v))
  low <- mu-sem
  mu[low<5] <- 0
  return(mu)
}

# Return SEM
standardError <- function(v){
  mu <- rowMeans(v)
  sem <- apply(v,1,sd)/sqrt(ncol(v))
  return(sem)
}

# Catenate data from different videosets
data <- read.csv("/path/data_avg_movieclips/average_ratings_videoset_1.csv")
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_2.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_3.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_4.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_5.csv"))
data <- rbind(data,read.csv("/path/data_avg_movieclips/average_ratings_videoset_6.csv"))

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
datapath <- "/path/data_csv_movieclips"
for(i in seq(1:length(data))){
  values <- c(values,data[,i])
  names <- c(names,rep(dims[i],234))
  mu <- c(mu,rep(mean(data[,i],na.rm=TRUE),234))
  
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
  alpha <- c(alpha,rep(mean(c(alphaV1$alpha,alphaV2$alpha,alphaV3$alpha,alphaV4$alpha,alphaV5$alpha,alphaV6$alpha)),234))
  
  # Occurrence rate
  occ <- occ_rate(v1)+occ_rate(v2)+occ_rate(v3)+occ_rate(v4)+occ_rate(v5)+occ_rate(v6)
  occurance <- c(occurance,rep(occ,234))
  
  # Mean values
  over_rate <- c(over(v1),over(v2),over(v3),over(v4),over(v5),over(v6))
  values_over <- c(values_over,over_rate)
  over_rate <- over_rate[over_rate>0]
  mu_over <- c(mu_over,rep(mean(over_rate,na.rm=TRUE),234))
  
  # SEM
  se <- c(se,rep(mean(c(standardError(v1),standardError(v2),standardError(v3),standardError(v4),standardError(v5),standardError(v6))),234))
  
  # Create a combined design matrix for ICC (absolute) calculation over all videosets
  m <- matrix(NA,234,sum(length(v1),length(v2),length(v3),length(v4),length(v5),length(v6)))
  m[1:39,1:length(v1)] <- as.matrix(v1)
  m[40:78,(length(v1)+1):sum(length(v1),length(v2))] <- as.matrix(v2)
  m[79:117,(sum(length(v1),length(v2))+1):sum(length(v1),length(v2),length(v3))] <- as.matrix(v3)
  m[118:156,(sum(length(v1),length(v2),length(v3))+1):sum(length(v1),length(v2),length(v3),length(v4))] <- as.matrix(v4)
  m[157:195,(sum(length(v1),length(v2),length(v3),length(v4))+1):sum(length(v1),length(v2),length(v3),length(v4),length(v5))] <- as.matrix(v5)
  m[196:234,(sum(length(v1),length(v2),length(v3),length(v4),length(v5))+1):sum(length(v1),length(v2),length(v3),length(v4),length(v5),length(v6))] <- as.matrix(v6)
  
  # Calculate ICCs
  res <- ICC(m,alpha=.05)
  
  # Single measures ICCs
  icc2 <- c(icc2,rep(res$results$ICC[2],234))
  icc2_lb <-c(icc2_lb,rep(res$results$`lower bound`[2],234))
  icc2_ub <-c(icc2_ub,rep(res$results$`upper bound`[2],234))
  
}

# Save a ICC data frame with raw rating values for plotting
df_icc2 <- data.frame(values,values_over,occurance,mu,se,mu_over,icc2,icc2_lb,icc2_ub,alpha)
df_icc2$occurance_d <- occurance/234
df_icc2$names <- factor(names,levels = unique(names))
#write.csv(df_icc2,file = "/path/consistency_measures_for_plotting.csv")

# Save also summary values for each video clip
ind <- seq(from=1,to=31824,by=234)
tbl_icc2 <- df_icc2[ind,seq(from=3,to=12)]
colnames(tbl_icc2) <- c("Occurrence","Mean values","SE","Mean values over 5","ICC","ICC_lb","ICC_ub","Alpha","Occurrence_d","Dimension")
#write.csv(tbl_icc2,file = "/path/consistency_measures.csv")