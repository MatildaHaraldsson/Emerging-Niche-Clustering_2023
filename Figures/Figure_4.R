#######################################################################################
# Script for making Figure 4 in Haraldsson and Thébault 2023, Ecology Letters
#######################################################################################

rm(list=ls())

# libraries
library(ggplot2)
library(viridis)

# load necessary data
(load("long_clust_comp_shuffle.RData"))   # nr of clusters - results from successful simulation and clustering analysis for model with competition only: nSim (simulation ID), time (time step), N.khat (number of clusters for competing community "N"), N.maxgap (clustering statistic for N), N.score (clustering statistic for N), N.pvalue (P.value for N)
(load("long_clust_corr_shuffle.RData"))   # nr of clusters - results from successful simulation and clustering analysis for correlated predator-prey model: nSim, time, N.khat (number of clusters for competing prey community "N"), N.maxgap, N.score, N.pvalue, P.khat (number of clusters for predator community "P"), P.maxgap, P.score, P.pvalue
(load("long_clust_uncorr_shuffle.RData")) # nr of clusters - results from successful simulation and clustering analysis for uncorrelated predator-prey model: nSim, time, N.khat (number of clusters for competing prey community on competition niche axis "N"), N.maxgap, N.score, N.pvalue, P.khat, P.maxgap, P.score, P.pvalue, M.khat (number of clusters in competing prey community on predator interaction niche axis "M"), M.maxgap, M.score, M.pvalue

(load("long_abu_comp.RData"))   # richness competing community - results from successful simulation for model with competition only: nSim, ts, nrSp 
(load("long_abu_corr_N.RData")) # richness competing prey community - results from successful simulation for correlated predator-prey model: nSim (simulation ID), ts (time step), nrSp (number of species survived at that timestep)
(load("long_abu_corr_P.RData")) # richness predator community - results from successful simulation for correlated predator-prey model: nSim, ts, nrSp

(load("long_abu_uncorr_N.RData")) # richness competing prey community - results from successful simulation for uncorrelated predator-prey model: nSim, ts, nrSp
(load("long_abu_uncorr_P.RData")) # richness predator community - results from successful simulation for uncorrelated predator-prey model: nSim, ts, nrSp

#################################################################
### a) Correlated case, and b) Uncorrelated case
#################################################################

# Organize data
#################################################################

# pick out data for each community, excluding NA´s
long_clust_corr_N <- long_clust_corr[-which(is.na(long_clust_corr$N.pvalue)),]
long_clust_corr_P <- long_clust_corr[-which(is.na(long_clust_corr$P.pvalue)),]
long_clust_uncorr_N <- long_clust_uncorr[-which(is.na(long_clust_uncorr$N.pvalue)),]
long_clust_uncorr_M <- long_clust_uncorr[-which(is.na(long_clust_uncorr$M.pvalue)),]
long_clust_uncorr_P <- long_clust_uncorr[-which(is.na(long_clust_uncorr$P.pvalue)),]

# Mean number of clusters (+/- SE) with time
#################################################################

# organize data
varName <- ".khat" # cluster variable
data_allNames <- c("long_clust_comp","long_clust_corr_N","long_clust_corr_P",
                   "long_clust_uncorr_N","long_clust_uncorr_M","long_clust_uncorr_P")

timeSteps <- sort(unique(long_clust_comp$time))
meanDat1 <- meanDat2 <- meanDat3 <- meanDat4 <- meanDat5 <- meanDat6 <- c()
for(i in 1:length(timeSteps)) {
  TheTimeStep <- timeSteps[i]
  for(j in 1:6){ # calculate statistics for each dataset 1) competing community, 2) prey community (correlated and uncorrelated), 3) predator community  (correlated and uncorrelated)
    data_all <- get(data_allNames[j])
    dat_select <- data_all[data_all$time == TheTimeStep,]
    
    N.name <- paste0("N",varName)
    P.name <- paste0("P",varName)
    M.name <- paste0("M",varName)
    
    theMean <- mean(dat_select[,N.name])
    theSE <- sd(dat_select[,N.name])/sqrt(length(dat_select[,N.name]))
    theMin <- theMean-theSE
    theMax <- theMean+theSE
    if(j == 1) meanDat1 <- rbind(meanDat1,c(theMean,theSE,theMin,theMax))
    if(j == 2) meanDat2 <- rbind(meanDat2,c(theMean,theSE,theMin,theMax))
    if(j == 3) {
      theMean <- mean(dat_select[,P.name])
      theSE <- sd(dat_select[,P.name]) /sqrt(length(dat_select[,P.name])) 
      theMin <- theMean-theSE
      theMax <- theMean+theSE
      meanDat3 <- rbind(meanDat3,c(theMean,theSE,theMin,theMax))
    } 
    if(j == 4) meanDat4 <- rbind(meanDat4,c(theMean,theSE,theMin,theMax))
    if(j == 5) {
      theMean <- mean(dat_select[,M.name])
      theSE <- sd(dat_select[,M.name]) /sqrt(length(dat_select[,M.name])) 
      theMin <- theMean-theSE
      theMax <- theMean+theSE
      meanDat5 <- rbind(meanDat5,c(theMean,theSE,theMin,theMax))
    }   
    if(j == 6) {
      theMean <- mean(dat_select[,P.name])
      theSE <- sd(dat_select[,P.name]) /sqrt(length(dat_select[,P.name])) 
      theMin <- theMean-theSE
      theMax <- theMean+theSE
      meanDat6 <- rbind(meanDat6,c(theMean,theSE,theMin,theMax))
    } # end if statement    
  } # end j loop
} # end i loop

meanDatMat <- cbind(meanDat1,meanDat2,meanDat3)
colnames(meanDatMat) <- c("mean1","se1","lower1","upper1","mean2","se2","lower2","upper2","mean3","se3","lower3","upper3")
meanDatMatFrame <- data.frame(time=timeSteps,meanDatMat)

meanDatMatComp <- cbind(meanDat1,meanDat4,meanDat5,meanDat6)
colnames(meanDatMatComp) <- c("mean1","se1","lower1","upper1","mean4","se4","lower4","upper4","mean5","se5","lower5","upper5","mean6","se6","lower6","upper6")
meanDatMatFrame2NicheComp <- data.frame(time=timeSteps,meanDatMatComp)
rm(meanDat1,meanDat2,meanDat3,meanDat4,meanDat5,meanDat6)

# make plotstructure - correlated case
data_all <- as.data.frame(meanDatMatFrame)
plotVar <- "Clusters"

tClust <- ggplot(data_all) + 
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 16)) +
  geom_line(aes(x=rep(500,59),y=seq(0,16,length.out=59)),linetype="dashed",color="darkgray") +
  geom_ribbon(aes(x=time,ymin=lower1,ymax=upper1,fill=viridis(3,begin=0,end=1,option="D")[3]),alpha=0.7) +
  geom_point(aes(x=time,y=mean1),size=1,color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_line(aes(x=time,y=mean1),color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_ribbon(aes(x=time,ymin=lower2,ymax=upper2,fill=viridis(3,begin=0,end=1,option="D")[2]),alpha=0.7) +
  geom_point(aes(x=time,y=mean2),size=1,color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_line(aes(x=time,y=mean2),color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_ribbon(aes(x=time,ymin=lower3,ymax=upper3,fill=viridis(3,begin=0,end=1,option="D")[1]),alpha=0.7) +
  geom_point(aes(x=time,y=mean3),size=1,color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  geom_line(aes(x=time,y=mean3),color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  scale_fill_manual(name="",    # NOTE - Order for colors different than above
                    values=c(viridis(3,begin=0,end=1,option="D")[2],viridis(3,begin=0,end=1,option="D")[1],viridis(3,begin=0,end=1,option="D")[3]),
                    labels=c("Prey","Predator","Competing community")) +
  scale_y_continuous(name=paste(plotVar), breaks=seq(0,16,length.out=3)) +
  scale_x_continuous(name="",breaks=seq(0,5000,length.out=3)) +
  theme_minimal(base_size=10, base_line_size=0.1) +
  theme(legend.position="bottom",legend.box="vertical",legend.text = element_text(size=8), 
        legend.title = element_text(size=8),legend.key.size = unit(0.20, 'cm'),
        plot.margin = unit(c(0,0.1,0,0.3), "cm"))



# make plotstructure - uncorrelated case
data_all <- as.data.frame(meanDatMatFrame2NicheComp)
plotVar <- "Clusters"

tClust2Niche <- ggplot(data_all) + 
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 16)) +
  geom_line(aes(x=rep(500,59),y=seq(0,16,length.out=59)),linetype="dashed",color="darkgray") +
  geom_ribbon(aes(x=time,ymin=lower1,ymax=upper1,fill=viridis(6,begin=0,end=1,option="D")[6]),alpha=0.7) +
  geom_point(aes(x=time,y=mean1),size=1,color=c(viridis(6,begin=0,end=1,option="D")[6])) +
  geom_line(aes(x=time,y=mean1),color=c(viridis(6,begin=0,end=1,option="D")[6])) +
  geom_ribbon(aes(x=time,ymin=lower4,ymax=upper4,fill=viridis(6,begin=0,end=1,option="D")[5]),alpha=0.7) +
  geom_point(aes(x=time,y=mean4),size=1,color=viridis(6,begin=0,end=1,option="D")[5]) +
  geom_line(aes(x=time,y=mean4),color=viridis(6,begin=0,end=1,option="D")[5]) +
  geom_ribbon(aes(x=time,ymin=lower5,ymax=upper5,fill=viridis(6,begin=0,end=1,option="D")[3]),alpha=0.7) +
  geom_point(aes(x=time,y=mean5),size=1,color=c(viridis(6,begin=0,end=1,option="D")[3])) +
  geom_line(aes(x=time,y=mean5),color=c(viridis(6,begin=0,end=1,option="D")[3])) +
  geom_ribbon(aes(x=time,ymin=lower6,ymax=upper6,fill=viridis(6,begin=0,end=1,option="D")[2]),alpha=0.7) +
  geom_point(aes(x=time,y=mean6),size=1,color=c(viridis(6,begin=0,end=1,option="D")[2])) +
  geom_line(aes(x=time,y=mean6),color=c(viridis(6,begin=0,end=1,option="D")[2])) +
  scale_fill_manual(name="",    # NOTE - Order for colors different than above
                    values=c(viridis(6,begin=0,end=1,option="D")[3],
                             viridis(6,begin=0,end=1,option="D")[2],
                             viridis(6,begin=0,end=1,option="D")[5],
                             viridis(6,begin=0,end=1,option="D")[6]),
                    labels=c(expression(paste("Prey ", italic(nu^{N}))),
                             "Predator",
                             expression(paste("Prey ", italic(mu^{N}))),
                             "Competing community") ) +
  scale_y_continuous(name=paste(plotVar), breaks=seq(0,16,length.out=3)) +
  scale_x_continuous(name="",breaks=seq(0,5000,length.out=3)) +
  theme_minimal(base_size=10, base_line_size=0.1) +
  theme(legend.position="bottom",legend.box="vertical",legend.text = element_text(size=8), 
        legend.title = element_text(size=8),legend.key.size = unit(0.20, 'cm'),
        plot.margin = unit(c(0,0.1,0,0.3), "cm"))

# Mean richness (+/- SE) with time
#################################################################
# Correlated case
varName <- "nrSp"
data_allNames <- c("long_abu_comp","long_abu_corr_N","long_abu_corr_P")

meanDat1 <- meanDat2 <-meanDat3 <- c()
for(i in 1:length(timeSteps)) {
  TheTimeStep <- timeSteps[i]
  for(j in 1:3){ # calculate for each dataset: 1) competing community, 2) prey community and 3) predator community for correlated case
    data_all <- get(data_allNames[j])
    dat_select <- data_all[data_all$ts == TheTimeStep,]
    theMean <- mean(na.omit(dat_select[,varName]))
    theSE <- sd(na.omit(dat_select[,varName])) /sqrt(length(na.omit(dat_select[,varName])))
    theMin <- theMean-theSE
    theMax <- theMean+theSE
    if(j == 1) meanDat1 <- rbind(meanDat1,c(theMean,theSE,theMin,theMax))
    if(j == 2) meanDat2 <- rbind(meanDat2,c(theMean,theSE,theMin,theMax))
    if(j == 3) meanDat3 <- rbind(meanDat3,c(theMean,theSE,theMin,theMax))
  }
}

meanDatVec <- cbind(meanDat1,meanDat2,meanDat3)
colnames(meanDatVec) <- c("mean1","se1","lower1","upper1","mean2","se2","lower2","upper2","mean3","se3","lower3","upper3")
meanDatFrame <- data.frame(time=timeSteps,meanDatVec)
rm(meanDat1,meanDat2,meanDat3)

# make plotstructure - uncorrelated case
data_all <- as.data.frame(meanDatFrame)
plotVar <- "Richness"
main_text <- c(paste(plotVar,"with time"))

tRich <- ggplot(data_all) + 
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 200)) +
  geom_line(aes(x=rep(500,dim(data_all)[1]),y=seq(0,200,length.out=dim(data_all)[1])),linetype="dashed",color="darkgray") +
  geom_ribbon(aes(x=time,ymin=lower1,ymax=upper1,fill=viridis(3,begin=0,end=1,option="D")[3]),alpha=0.7) +
  geom_point(aes(x=time,y=mean1),size=1,color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_line(aes(x=time,y=mean1),color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_ribbon(aes(x=time,ymin=lower2,ymax=upper2,fill=viridis(3,begin=0,end=1,option="D")[2]),alpha=0.7) +
  geom_point(aes(x=time,y=mean2),size=1,color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_line(aes(x=time,y=mean2),color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_ribbon(aes(x=time,ymin=lower3,ymax=upper3,fill=viridis(3,begin=0,end=1,option="D")[1]),alpha=0.7) +
  geom_point(aes(x=time,y=mean3),size=1,color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  geom_line(aes(x=time,y=mean3),color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  scale_fill_manual(name="Community",    # NOTE - Order for colors different than above
                    values=c(viridis(3,begin=0,end=1,option="D")[2],viridis(3,begin=0,end=1,option="D")[1],viridis(3,begin=0,end=1,option="D")[3]),
                    labels=c("Prey","Pred","Competing")) +
  scale_y_continuous(name=paste(plotVar), breaks=seq(0,200,length.out=3)) +
  scale_x_continuous(name="",breaks=seq(0,5000,length.out=3)) +
  guides(fill="none") +
  theme_minimal(base_size=10, base_line_size=0.1)

# Uncorrelated case
varName <- "nrSp"
data_allNames <- c("long_abu_comp","long_abu_uncorr_N","long_abu_uncorr_P")

meanDat1 <- meanDat2 <-meanDat3 <- c()
for(i in 1:length(timeSteps)) {
  TheTimeStep <- timeSteps[i]
  for(j in 1:3){ # calculate for each dataset: 1) competing community, 2) prey community and 3) predator community for uncorrelated case
    data_all <- get(data_allNames[j])
    dat_select <- data_all[data_all$ts == TheTimeStep,]
    theMean <- mean(na.omit(dat_select[,varName]))
    theSE <- sd(na.omit(dat_select[,varName])) /sqrt(length(na.omit(dat_select[,varName])))
    theMin <- theMean-theSE
    theMax <- theMean+theSE
    if(j == 1) meanDat1 <- rbind(meanDat1,c(theMean,theSE,theMin,theMax))
    if(j == 2) meanDat2 <- rbind(meanDat2,c(theMean,theSE,theMin,theMax))
    if(j == 3) meanDat3 <- rbind(meanDat3,c(theMean,theSE,theMin,theMax))
  }
}

meanDatVec <- cbind(meanDat1,meanDat2,meanDat3)
colnames(meanDatVec) <- c("mean1","se1","lower1","upper1","mean2","se2","lower2","upper2","mean3","se3","lower3","upper3")
meanDatFrame <- data.frame(time=timeSteps,meanDatVec)
rm(meanDat1,meanDat2,meanDat3)


# make plotstructure - uncorrelated case
data_all <- as.data.frame(meanDatFrame)
plotVar <- "Richness"
main_text <- c(paste(plotVar,"with time"))

tRich2Niche <- ggplot(data_all) + 
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 200)) +
  geom_line(aes(x=rep(500,dim(data_all)[1]),y=seq(0,200,length.out=dim(data_all)[1])),linetype="dashed",color="darkgray") +
  geom_ribbon(aes(x=time,ymin=lower1,ymax=upper1,fill=viridis(3,begin=0,end=1,option="D")[3]),alpha=0.7) +
  geom_point(aes(x=time,y=mean1),size=1,color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_line(aes(x=time,y=mean1),color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_ribbon(aes(x=time,ymin=lower2,ymax=upper2,fill=viridis(3,begin=0,end=1,option="D")[2]),alpha=0.7) +
  geom_point(aes(x=time,y=mean2),size=1,color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_line(aes(x=time,y=mean2),color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_ribbon(aes(x=time,ymin=lower3,ymax=upper3,fill=viridis(3,begin=0,end=1,option="D")[1]),alpha=0.7) +
  geom_point(aes(x=time,y=mean3),size=1,color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  geom_line(aes(x=time,y=mean3),color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  scale_fill_manual(name="Community",    # NOTE - Order for colors different than above
                    values=c(viridis(3,begin=0,end=1,option="D")[2],viridis(3,begin=0,end=1,option="D")[1],viridis(3,begin=0,end=1,option="D")[3]),
                    labels=c("Prey","Pred","Competing")) +
  scale_y_continuous(name=paste(plotVar), breaks=seq(0,200,length.out=3)) +
  scale_x_continuous(name="",breaks=seq(0,5000,length.out=3)) +
  guides(fill="none") +
  theme_minimal(base_size=10, base_line_size=0.1)

# Mean richness per cluster (+/- SE) with time
#################################################################
# Correlated case

data_allNamesClust <- c("long_clust_comp","long_clust_corr_N","long_clust_corr_P")
data_allNamesAbu <- c("long_abu_comp","long_abu_corr_N","long_abu_corr_P")

varName <- "nrSp"
meanDat1 <- meanDat2 <-meanDat3 <- meanDat01 <- meanDat02 <-meanDat03 <- c()

for(i in 1:length(timeSteps)) {
  TheTimeStep <- timeSteps[i]
  
  for(j in 1:3){ # calculate for each dataset
    data_all <- get(data_allNamesClust[j])
    data_all_sim <- get(data_allNamesAbu[j])
    dat_select <- data_all[data_all$time == TheTimeStep,]
    dat_select_sim <- data_all_sim[data_all_sim$ts == TheTimeStep,]
    RichClust <- c()
    for(k in 1:length(dat_select$nSim)) {
      nSim <- dat_select[k,"nSim"]
      if(j != 3) meanRich <- dat_select_sim[which(dat_select_sim[,"nSim"] == nSim),"nrSp"] / dat_select[k,"N.khat"]
      if(j == 3) meanRich <- dat_select_sim[which(dat_select_sim[,"nSim"] == nSim),"nrSp"] / dat_select[k,"P.khat"]
      RichClust <- c(RichClust,meanRich)
    }
    # calculate mean for all data and data > 0
    theMean0 <- mean(RichClust)
    theMean <- mean(RichClust[RichClust > 0])
    theSE0 <- sd(RichClust)/sqrt(length(RichClust))
    theSE <- sd(RichClust[RichClust > 0])/sqrt(length(RichClust[RichClust > 0]))
    theMin0 <- theMean0-theSE0
    theMax0 <- theMean0+theSE0
    theMin <- theMean-theSE
    theMax <- theMean+theSE
    if(j == 1) {
      meanDat1 <- rbind(meanDat1,c(theMean,theSE,theMin,theMax)) 
      meanDat01 <- rbind(meanDat01,c(theMean0,theSE0,theMin0,theMax0))
    }
    if(j == 2) { 
      meanDat2 <- rbind(meanDat2,c(theMean,theSE,theMin,theMax)) 
      meanDat02 <- rbind(meanDat02,c(theMean0,theSE0,theMin0,theMax0))
    }
    if(j == 3) {
      meanDat3 <- rbind(meanDat3,c(theMean,theSE,theMin,theMax)) 
      meanDat03 <- rbind(meanDat03,c(theMean0,theSE0,theMin0,theMax0))
    }
  }
}

meanDatMat <- cbind(meanDat1,meanDat2,meanDat3)
#meanDatMat <- cbind(meanDat01,meanDat02,meanDat03)
colnames(meanDatMat) <- c("mean1","se1","lower1","upper1","mean2","se2","lower2","upper2","mean3","se3","lower3","upper3")
meanDatMatFrame <- data.frame(time=timeSteps,meanDatMat)
rm(meanDat1,meanDat2,meanDat3)

data_all <- as.data.frame(meanDatMatFrame)
plotVar <-    str_wrap("Richness within clusters", width = 15)

tRichClust <- ggplot(data_all) + 
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 100)) +
  geom_line(aes(x=rep(500,59),y=seq(0,100,length.out=59)),linetype="dashed",color="darkgray") +
  geom_ribbon(aes(x=time,ymin=lower1,ymax=upper1,fill=viridis(3,begin=0,end=1,option="D")[3]),alpha=0.7) +
  geom_point(aes(x=time,y=mean1),size=1,color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_line(aes(x=time,y=mean1),color=c(viridis(3,begin=0,end=1,option="D")[3])) +
  geom_ribbon(aes(x=time,ymin=lower2,ymax=upper2,fill=viridis(3,begin=0,end=1,option="D")[2]),alpha=0.7) +
  geom_point(aes(x=time,y=mean2),size=1,color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_line(aes(x=time,y=mean2),color=c(viridis(3,begin=0,end=1,option="D")[2])) +
  geom_ribbon(aes(x=time,ymin=lower3,ymax=upper3,fill=viridis(3,begin=0,end=1,option="D")[1]),alpha=0.7) +
  geom_point(aes(x=time,y=mean3),size=1,color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  geom_line(aes(x=time,y=mean3),color=c(viridis(3,begin=0,end=1,option="D")[1])) +
  scale_fill_manual(name="Community",    # NOTE - Order for colors different than above
                    values=c(viridis(3,begin=0,end=1,option="D")[2],viridis(3,begin=0,end=1,option="D")[1],viridis(3,begin=0,end=1,option="D")[3]),
                    labels=c("Prey","Pred","Competing")) +
  scale_y_continuous(name=paste(plotVar), breaks=seq(0,100,length.out=3)) +
  scale_x_continuous(name="",breaks=seq(0,5000,length.out=3)) +
  guides(fill="none") +
  theme_minimal(base_size=10, base_line_size=0.1)

# Uncorrelated case

data_allNamesClust <- c("long_clust_comp","long_clust_uncorr_N","long_clust_uncorr_M","long_clust_uncorr_P")
data_allNamesAbu <- c("long_abu_comp","long_abu_uncorr_N","long_abu_uncorr_N","long_abu_uncorr_P")

varName <- "nrSp"

meanDat1 <- meanDat2 <- meanDat3 <- meanDat4 <- meanDat01 <- meanDat02 <- meanDat03 <- meanDat04 <- c()

for(i in 1:length(timeSteps)) {
  TheTimeStep <- timeSteps[i]
  
  for(j in 1:4){ # calculate for each dataset
    data_all <- get(data_allNamesClust[j])
    data_all_sim <- get(data_allNamesAbu[j])
    dat_select <- data_all[data_all$time == TheTimeStep,]
    dat_select_sim <- data_all_sim[data_all_sim$ts == TheTimeStep,]
    RichClust <- c()
    for(k in 1:length(dat_select$nSim)) {
      nSim <- dat_select[k,"nSim"]
      if(j == 1 | j == 2) meanRich <- dat_select_sim[which(dat_select_sim[,"nSim"] == nSim),"nrSp"] / dat_select[k,"N.khat"]
      if(j == 3) meanRich <- dat_select_sim[which(dat_select_sim[,"nSim"] == nSim),"nrSp"] / dat_select[k,"M.khat"]
      if(j == 4) meanRich <- dat_select_sim[which(dat_select_sim[,"nSim"] == nSim),"nrSp"] / dat_select[k,"P.khat"]
      RichClust <- c(RichClust,meanRich)
    }
    # calculate mean for all data and data > 0
    theMean0 <- mean(RichClust)
    theMean <- mean(RichClust[RichClust > 0])
    theSE0 <- sd(RichClust)/sqrt(length(RichClust))
    theSE <- sd(RichClust[RichClust > 0])/sqrt(length(RichClust[RichClust > 0]))
    theMin0 <- theMean0-theSE0
    theMax0 <- theMean0+theSE0
    theMin <- theMean-theSE
    theMax <- theMean+theSE
    if(j == 1) {
      meanDat1 <- rbind(meanDat1,c(theMean,theSE,theMin,theMax)) 
      meanDat01 <- rbind(meanDat01,c(theMean0,theSE0,theMin0,theMax0))
    }
    if(j == 2) { 
      meanDat2 <- rbind(meanDat2,c(theMean,theSE,theMin,theMax)) 
      meanDat02 <- rbind(meanDat02,c(theMean0,theSE0,theMin0,theMax0))
    }
    if(j == 3) {
      meanDat3 <- rbind(meanDat3,c(theMean,theSE,theMin,theMax)) 
      meanDat03 <- rbind(meanDat03,c(theMean0,theSE0,theMin0,theMax0))
    }
  if(j == 4) {
    meanDat4 <- rbind(meanDat4,c(theMean,theSE,theMin,theMax)) 
    meanDat04 <- rbind(meanDat04,c(theMean0,theSE0,theMin0,theMax0))
    }
  }
} 

meanDatMat <- cbind(meanDat1,meanDat2,meanDat3,meanDat4)
colnames(meanDatMat) <- c("mean1","se1","lower1","upper1","mean2","se2","lower2","upper2","mean3","se3","lower3","upper3","mean4","se4","lower4","upper4")
meanDatMatFrame <- data.frame(time=timeSteps,meanDatMat)
rm(meanDat1,meanDat2,meanDat3)

data_all <- as.data.frame(meanDatMatFrame)
plotVar <-    str_wrap("Richness within clusters", width = 15)

tRichClust2Niche <- ggplot(data_all) + 
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 100)) +
  geom_line(aes(x=rep(500,59),y=seq(0,100,length.out=59)),linetype="dashed",color="darkgray") +
  geom_point(aes(x=time,y=mean1),size=1,color=c(viridis(6,begin=0,end=1,option="D")[6])) +
  geom_line(aes(x=time,y=mean1),color=c(viridis(6,begin=0,end=1,option="D")[6])) +
  geom_ribbon(aes(x=time,ymin=lower2,ymax=upper2,fill=viridis(6,begin=0,end=1,option="D")[5]),alpha=0.7) +
  geom_point(aes(x=time,y=mean2),size=1,color=c(viridis(6,begin=0,end=1,option="D")[5])) +
  geom_line(aes(x=time,y=mean2),color=c(viridis(6,begin=0,end=1,option="D")[5])) +
  geom_ribbon(aes(x=time,ymin=lower3,ymax=upper3,fill=viridis(6,begin=0,end=1,option="D")[3]),alpha=0.7) +
  geom_point(aes(x=time,y=mean3),size=1,color=c(viridis(6,begin=0,end=1,option="D")[3])) +
  geom_line(aes(x=time,y=mean3),color=c(viridis(6,begin=0,end=1,option="D")[3])) +
  geom_ribbon(aes(x=time,ymin=lower4,ymax=upper4,fill=viridis(6,begin=0,end=1,option="D")[2]),alpha=0.7) +
  geom_point(aes(x=time,y=mean4),size=1,color=c(viridis(6,begin=0,end=1,option="D")[2])) +
  geom_line(aes(x=time,y=mean4),color=c(viridis(6,begin=0,end=1,option="D")[2])) +
  scale_fill_manual(name="Community",    # NOTE - Order for colors different than above
                    values=c(viridis(6,begin=0,end=1,option="D")[3],
                             viridis(6,begin=0,end=1,option="D")[2],
                             viridis(6,begin=0,end=1,option="D")[5],
                             viridis(6,begin=0,end=1,option="D")[6]),
                    labels=c("Prey","Pred","Competing","Competing2")) +
  scale_y_continuous(name=paste(plotVar), breaks=seq(0,100,length.out=3)) +
  scale_x_continuous(name="",breaks=seq(0,5000,length.out=3)) +
  guides(fill="none") +
  theme_minimal(base_size=10, base_line_size=0.1)

#################################################################
### make the figure
#################################################################

# save in tiff format
#####################
tiff("Fig4.tiff", units="in", width=(5), height=3.9, res=300)

layout <- "
#ABBBBCCCCDDDD
#EFFFFGGGGHHHH
#EFFFFGGGGHHHH
#EFFFFGGGGHHHH
#EFFFFGGGGHHHH
#EMMMMNNNNOOOO
#EIIIIJJJJKKKK
#EIIIIJJJJKKKK
#EIIIIJJJJKKKK
#EIIIIJJJJKKKK
"
plot_spacer() + plot_spacer() +plot_spacer() + plot_spacer() +
  plot_spacer() + tRich + tClust + tRichClust +
  tRich2Niche + tClust2Niche + tRichClust2Niche +
  plot_spacer() +plot_spacer() + plot_spacer() +
  plot_layout(design = layout)

grid.text("Generations",x=unit(0.53, "npc"),y=unit(0.15,"npc"),gp=gpar(fontsize=9, col="black")) 
grid.text("Generations",x=unit(0.53, "npc"),y=unit(0.63,"npc"),gp=gpar(fontsize=9, col="black")) 
grid.text("a) Correlated case",x=unit(0.25, "npc"),y=unit(0.95,"npc"),gp=gpar(fontsize=10, col="black")) 
grid.text("b) Uncorrelated case",x=unit(0.27, "npc"),y=unit(0.47,"npc"),gp=gpar(fontsize=10, col="black")) 

dev.off()

# save in pdf format
#####################
dir_plot <- paste0(getwd(),"/Fig4.pdf") # save the plot to the current work directory
pdf(file = dir_plot,   # filename and chosen directory
    width = 4.5,       # the width of the plot in inches
    height = 3.6      # the height of the plot in inches
)

layout <- "
#ABBBBCCCCDDDD
#EFFFFGGGGHHHH
#EFFFFGGGGHHHH
#EFFFFGGGGHHHH
#EFFFFGGGGHHHH
#EMMMMNNNNOOOO
#EIIIIJJJJKKKK
#EIIIIJJJJKKKK
#EIIIIJJJJKKKK
#EIIIIJJJJKKKK
"
plot_spacer() + plot_spacer() +plot_spacer() + plot_spacer() +
  plot_spacer() + tRich + tClust + tRichClust +
  tRich2Niche + tClust2Niche + tRichClust2Niche +
  plot_spacer() +plot_spacer() + plot_spacer() +
  plot_layout(design = layout)

grid.text("Generations",x=unit(0.53, "npc"),y=unit(0.15,"npc"),gp=gpar(fontsize=9, col="black")) 
grid.text("Generations",x=unit(0.53, "npc"),y=unit(0.63,"npc"),gp=gpar(fontsize=9, col="black")) 
grid.text("a) Correlated case",x=unit(0.25, "npc"),y=unit(0.95,"npc"),gp=gpar(fontsize=10, col="black")) 
grid.text("b) Uncorrelated case",x=unit(0.27, "npc"),y=unit(0.47,"npc"),gp=gpar(fontsize=10, col="black")) 

dev.off()