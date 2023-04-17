#######################################################################################
# Script for making Figure 2 and 3 in Haraldsson and Thébault 2023, Ecology Letters
#######################################################################################

rm(list=ls())

# libraries
library(ggplot2)
library(viridis)
#library(tidyverse)
#library(scales)
#library(patchwork)
#library(gridExtra)
#library(grid)
#library(RGraphics)
source("Library_Emerging_niche_clustering.R")

# load simulation results
(load("datCorr.R"))
(load("datUncorr.R"))

set.seed(1)
#################################################################
### Figure 2 Correlated case 
#################################################################

# Correlation of niche axes
#################################################################

data_all <- as.data.frame(cbind(N=rep((runif(200,0,1))),M=(runif(200,0,1))))

corr_case <- ggplot(data_all) + 
  geom_point(aes(x=N,y=N),size=1,color="#21908CFF") +
  scale_y_continuous(breaks=seq(0,1,length.out=3)) +
  scale_x_continuous(name=expression(paste("Prey niche axis ", italic(mu^{N}) )),breaks=seq(0,1,length.out=3)) +
  theme(axis.text=element_text(size=7),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))


#################################################################
### a) Prey community
#################################################################

# Competition and predation kernels
#################################################################
mu_N <- seq(0,1,1/100) # niche axis for generating alpha values
a_N <- alpha_N(mu=mu_N,sigma_sp=0.10,nr_sp=length(mu_N),Lmax=1,P=4,Sigma_term=1) # generate alpha values for plotting competition kernel
colnames(a_N) <- round(mu_N,digits=3)
rownames(a_N) <- round(mu_N,digits=3)

dat <- as.data.frame(cbind(a=a_N[,"0.4"],mu=mu_N)) # target species at mu = 0.4 in the competition matrix
datPointN <- as.data.frame(cbind(a=rep(0,length(sort(runif(20,0,1)))),mu=sort(runif(20,0,1))))
datPointP <- as.data.frame(cbind(a=rep(0,length(sort(runif(20,0,1)))),mu=sort(runif(20,0,1))))

# prey niche axis
preyniche1 <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_vline(xintercept=0.4, linetype="dotdash",color="darkgray") +
  geom_line(data=dat, colour="darkgray") +
  geom_point(data=datPointN, shape=21, size=1, color="#21908CFF") +
  geom_point(aes(x=0.4,y=0), size=3, color="#21908CFF") +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(expression(paste("Prey niche axis ", italic(mu^{N})))) +
  geom_segment(aes(x=0.4,y=0.5,xend=0.49,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.49,y=0.5,xend=0.4,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.4,y=-0.4,xend=0.4,yend=0),arrow=arrow(angle=9,length=unit(0.10,"inches"),type="closed")) +
  annotate("text", x=0.62, y=0.7, label = expression(paste(italic(sigma^{N}))),size=4) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))

# predator niche axis
predniche1 <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_vline(xintercept=0.4, linetype="dotdash",color="darkgray") +
  geom_line(data=dat, colour="darkgray") +
  geom_point(data=datPointN, shape=21, size=1, color="#21908CFF") +  
  coord_cartesian(xlim = c(0,1),ylim = c(0,1.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(expression(paste("Prey niche axis ", italic(nu^{N})))) +
  geom_segment(aes(x=0.4,y=0.5,xend=0.49,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.49,y=0.5,xend=0.4,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  annotate("text", x=0.62, y=0.7, label = expression(paste(italic(sigma^{P}))),size=4) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))

# predator top axis
predUpniche1 <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_point(aes(x=0.4,y=0),size=3,color="#440154FF") +
  geom_point(data=datPointP, shape=21, size=1, color="#440154FF") +
  coord_cartesian(xlim = c(0,1),ylim = c(0,0.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(NULL) +
  geom_segment(aes(x=0.4,y=0,xend=0.4,yend=-0.15),arrow=arrow(angle=9,length=unit(0.10,"inches"),type="closed")) +
  annotate("text", x=0.55, y=0.14, label = expression(paste("Predator niche axis ",italic(mu^{P}))),size=3.5) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))

# Network - focus on prey community - in Figure 2 and 3
#################################################################
datA <- data.frame(x=c(0,3,6),y=c(2.95,2.95,2.95))
datB <- data.frame(x=c(0,3,6),y=c(1,1,1))

net1 <- ggplot(NULL,aes(x=x, y=y)) + 
  theme_void() +
  scale_fill_manual(values=c("darkgray","#21908CFF")) +
  geom_segment(aes(x=0,y=1,xend=0,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=1,xend=3,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=1,xend=6,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=3,y=1,xend=0,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=1,xend=3,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=1,xend=6,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=6,y=1,xend=0,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=1,xend=3,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=1,xend=6,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_point(data=datA,shape=21,aes(fill="datA"),size=8,show.legend = FALSE) +
  geom_point(data=datB,shape=21,aes(fill="datB"),size=8,show.legend = FALSE) +
  coord_cartesian(xlim = c(0,6), ylim = c(1,5), clip = "off")


#################################################################
### b) Predator community
#################################################################

# Competition and predation kernels
#################################################################

# predator niche axis
predniche2 <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_vline(xintercept=0.4, linetype="dotdash",color="darkgray") +
  geom_line(data=dat, colour="darkgray") +
  geom_point(data=datPointN, shape=21, size=1, color="#21908CFF") +  
  coord_cartesian(xlim = c(0,1),ylim = c(0,1.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(expression(paste("Prey niche axis ", italic(nu^{N})))) +
  geom_segment(aes(x=0.4,y=0.5,xend=0.49,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.49,y=0.5,xend=0.4,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.4,y=1,xend=0.4,yend=1.5),arrow=arrow(angle=9,length=unit(0.10,"inches"),type="closed")) +
  annotate("text", x=0.62, y=0.7, label = expression(paste(italic(sigma^{P}))),size=4) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))

# predator top axis
predUpniche2 <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_point(aes(x=0.4,y=0),size=3,color="#440154FF") +
  geom_point(data=datPointP, shape=21, size=1, color="#440154FF") +
  coord_cartesian(xlim = c(0,1),ylim = c(0,0.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(NULL) +
  annotate("text", x=0.55, y=0.14, label = expression(paste("Predator niche axis ",italic(mu^{P}))),size=3.5) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))



# Network - focus on predator community - in Figure 2 and 3
#################################################################
datA <- data.frame(x=c(0,3,6),y=c(2.95,2.95,2.95))
datB <- data.frame(x=c(0,3,6),y=c(1,1,1))

net2 <- ggplot(NULL,aes(x=x, y=y)) + 
  theme_void() +
  scale_fill_manual(values=c("#453781FF","darkgray")) +
  geom_segment(aes(x=0,y=1,xend=0,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=1,xend=3,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=1,xend=6,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=3,y=1,xend=0,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=1,xend=3,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=1,xend=6,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=6,y=1,xend=0,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=1,xend=3,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=1,xend=6,yend=2.47),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_point(data=datA,shape=21,aes(fill="datA"),size=8,show.legend = FALSE) +
  geom_point(data=datB,shape=21,aes(fill="datB"),size=8,show.legend = FALSE) +
  coord_cartesian(xlim = c(0,6), ylim = c(1,5), clip = "off")


# Density plots - prey and predator community (bottom of panels)
#################################################################
# plotting parameters
P.sign <- datCorr$P.pvalue <= 0.05
N.sign <- datCorr$N.pvalue <= 0.05
dat125 <- cbind(datCorr,N.sign,P.sign)
dat_code <- c(125)
K_range <- c(5,500)
com_range <- c("N","P")
com_name <- c("prey","predator")
marg1 <- marg2 <- marg3 <- marg4 <- 0.25

# Significant results %
#################################################################

for(i in NROW(dat_code)) { # for each dataset
  for(j in 1:NROW(K_range)) { # for each unique value of K
    for(k in 1:NROW(com_range)){ # for prey and predator community
      dat_all <- get(paste0("dat",dat_code[i]))
      x.var <- "sigma_spN"
      y.var <- "sigma_spP"
      z.var <- paste0(com_range[k],".sign")
      x.name <- "sigma N"
      y.name <- "sigma P"
      z.name <- "significant"
      dat_uniq_d <- as.data.frame(subset(dat_all,dat_all[,"K"] == K_range[j]))
    
      x <- unique(dat_uniq_d[,x.var])[order(unique(dat_uniq_d[,x.var]))]
      y <- unique(dat_uniq_d[,y.var])[order(unique(dat_uniq_d[,y.var]))]
      df <- expand.grid(x,y)
      z <- sapply(c(1:dim(df)[1]),function(l){
        var.1 <- df[l,1]
        var.2 <- df[l,2]
        z.crit <- (dat_uniq_d[,x.var] == var.1 & dat_uniq_d[,y.var] == var.2)
        z <- sum(na.omit(dat_uniq_d[z.crit,z.var]))/length((dat_uniq_d[z.crit,z.var]))
        return(z)
      })
      df$z <- z
      colnames(df) <- c("x","y","z")
      
      if(j == 1) assign(x=paste0("SI_sp",k,"_K",j),ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(0,1)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill="none"))
      
      if(j == 2) assign(x=paste0("SI_sp",k,"_K",j), ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(0,1),breaks=c(0,1)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.y=element_blank(),axis.text.x = element_text(vjust = -1),legend.position="bottom",legend.box="vertical",legend.text = element_text(size=7),  # different!
                                legend.title = element_text(size=0),legend.key.size = unit(0.45, 'cm'),
                                plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +  
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          
                          guides(fill = guide_colorbar(title = " ")))
    } 
  } 
} 

# Number of Clusters
#################################################################

for(i in NROW(dat_code)) { # for each dataset
  for(j in 1:NROW(K_range)) { # for each unique value of K
    for(k in 1:NROW(com_range)){ # for prey and predator community
      dat_all <- get(paste0("dat",dat_code[i]))
      dat_all <- dat_all[which(dat_all[,paste0(com_range[k],".pvalue")] <= 0.05),]
      x.var <- "sigma_spN"
      y.var <- "sigma_spP"
      z.var <- paste0(com_range[k],".khat")
      if(k == 1) x.name <- "" else x.name <- "Sigma N"
      y.name <- "Sigma P"
      z.name <- "Clusters"
      dat_uniq_d <- as.data.frame(subset(dat_all,dat_all[,"K"] == K_range[j])) 

      x <- unique(dat_uniq_d[,x.var])[order(unique(dat_uniq_d[,x.var]))]
      y <- unique(dat_uniq_d[,y.var])[order(unique(dat_uniq_d[,y.var]))]
      df <- expand.grid(x,y)
      z <- sapply(c(1:dim(df)[1]),function(l){
        var.1 <- df[l,1]
        var.2 <- df[l,2]
        z.crit <- (dat_uniq_d[,x.var] == var.1 & dat_uniq_d[,y.var] == var.2)
        z <- mean(na.omit(dat_uniq_d[z.crit,z.var]))
        return(z)
      })
      df$z <- z
      colnames(df) <- c("x","y","z")
    
      if(j == 1) assign(x=paste0("CL_sp",k,"_K",j),ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(1,17.5)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.x=element_blank(), plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill="none"))
      
      if(j == 2) assign(x=paste0("CL_sp",k,"_K",j), ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(1,17.5),breaks=c(1,17)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.y = element_text(hjust = 20), axis.text.x = element_text(vjust = -1),
                                legend.position="bottom",legend.box="vertical",legend.text = element_text(size=7), 
                                legend.title = element_text(size=0),legend.key.size = unit(0.45, 'cm'),
                                plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
                          annotate("text", x = -0.03, y = 0.36, angle = 90, label = expression(paste( italic(sigma^{P}))), size=4) +
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill = guide_colorbar(title = " ")))
    } 
  }
}


#################################################################
### Figure 3 Uncorrelated case
#################################################################

# Correlation of niche axes
#################################################################

uncorr_case <- ggplot(data_all) + 
  geom_point(aes(x=N,y=M),size=1,color="#73D055FF") +
  scale_y_continuous(breaks=seq(0,1,length.out=3)) +
  scale_x_continuous(name=expression(paste("Prey niche axis ", italic(mu^{N}) )),breaks=seq(0,1,length.out=3)) +
  theme(axis.text=element_text(size=7),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))

#################################################################
### a) Prey community
#################################################################

# Competition and predation kernels
#################################################################

datPointM <- as.data.frame(cbind(a=rep(0,length(sort(runif(20,0,1)))),mu=sort(runif(20,0,1))))

# prey niche axis
prey2niche <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_vline(xintercept=0.4, linetype="dotdash",color="darkgray") +
  geom_line(data=dat, colour="darkgray") +
  geom_point(data=datPointN, shape=21, size=1, color="#73D055FF") +
  geom_point(aes(x=0.4,y=0), size=3, color="#73D055FF") +
  coord_cartesian(xlim = c(0,1),ylim = c(0,1.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(expression(paste("Prey niche axis ", italic(mu^{N})))) +
  geom_segment(aes(x=0.4,y=0.5,xend=0.49,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.49,y=0.5,xend=0.4,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.4,y=-0.4,xend=0.4,yend=0),arrow=arrow(angle=9,length=unit(0.10,"inches"),type="closed")) +
  annotate("text", x=0.62, y=0.7, label = expression(paste(italic(sigma^{N}))),size=4) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))

# predator niche axis
pred2niche <- ggplot(NULL,aes(y=a,x=mu)) +
  geom_vline(xintercept=0.4, linetype="dotdash",color="darkgray") +
  geom_line(data=dat, colour="darkgray") +
  geom_point(data=datPointM, shape=21, size=1, color="#33638DFF") +  
  coord_cartesian(xlim = c(0,1),ylim = c(0,1.2), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(expression(paste("Prey niche axis ", italic(nu^{N})))) +
  geom_segment(aes(x=0.4,y=0.5,xend=0.49,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.49,y=0.5,xend=0.4,yend=0.5), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  annotate("text", x=0.62, y=0.7, label = expression(paste(italic(sigma^{P}))),size=4) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=10),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"))


# Density plots - prey and predator community (bottom of panels)
#################################################################

# plotting parameters
P.sign <- datUncorr$P.pvalue <= 0.05
N.sign <- datUncorr$N.pvalue <= 0.05
M.sign <- datUncorr$M.pvalue <= 0.05
dat225 <- cbind(datUncorr,N.sign,M.sign,P.sign)
dat_code <- c(225)
K_range <- c(5,500)
com_range <- c("N","P","M")
com_name <- c("prey","predator","prey with pred niche")

# Significant results %
#################################################################
for(i in NROW(dat_code)) { # for each dataset
  for(j in 1:NROW(K_range)) { # for each unique value of K
    for(k in 1:NROW(com_range)){ # for prey and predator community
      
      dat_all <- get(paste0("dat",dat_code[i]))
      x.var <- "sigma_spN"
      y.var <- "sigma_spP"
      z.var <- paste0(com_range[k],".sign")
      x.name <- "sigma N"
      y.name <- "sigma P"
      z.name <- "significant"
      dat_uniq_d <- as.data.frame(subset(dat_all,dat_all[,"K"] == K_range[j])) 

      x <- unique(dat_uniq_d[,x.var])[order(unique(dat_uniq_d[,x.var]))]
      y <- unique(dat_uniq_d[,y.var])[order(unique(dat_uniq_d[,y.var]))]
      df <- expand.grid(x,y)
      z <- sapply(c(1:dim(df)[1]),function(l){
        var.1 <- df[l,1]
        var.2 <- df[l,2]
        z.crit <- (dat_uniq_d[,x.var] == var.1 & dat_uniq_d[,y.var] == var.2)  
        z <- sum(na.omit(dat_uniq_d[z.crit,z.var]))/length((dat_uniq_d[z.crit,z.var]))
        return(z)
      })
      df$z <- z
      colnames(df) <- c("x","y","z")
      
      if(j == 1) assign(x=paste0("SI2_sp",k,"_K",j),ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(0,1)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill="none"))
      
      if(j == 2) assign(x=paste0("SI2_sp",k,"_K",j), ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(0,1),breaks=c(0,1)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.y=element_blank(),axis.text.x = element_text(vjust = -1),legend.position="bottom",legend.box="vertical",legend.text = element_text(size=7),  # different!
                                legend.title = element_text(size=0),legend.key.size = unit(0.45, 'cm'),
                                plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) + 
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill = guide_colorbar(title = " ")))
    }
  }
}

# Number of Clusters
#################################################################
for(i in NROW(dat_code)) { # for each dataset
  for(j in 1:NROW(K_range)) { # for each unique value of K
    for(k in 1:NROW(com_range)){ # for prey and predator community
      
      dat_all <- get(paste0("dat",dat_code[i]))
      dat_all <- dat_all[which(dat_all[,paste0(com_range[k],".pvalue")] <= 0.05),]
      x.var <- "sigma_spN"
      y.var <- "sigma_spP"
      z.var <- paste0(com_range[k],".khat")
      if(k == 1) x.name <- "" else x.name <- "Sigma N"
      y.name <- "Sigma P"
      z.name <- "Clusters"
      dat_uniq_d <- as.data.frame(subset(dat_all,dat_all[,"K"] == K_range[j]))
      
      x <- unique(dat_uniq_d[,x.var])[order(unique(dat_uniq_d[,x.var]))]
      y <- unique(dat_uniq_d[,y.var])[order(unique(dat_uniq_d[,y.var]))]
      df <- expand.grid(x,y)
      z <- sapply(c(1:dim(df)[1]),function(l){
        var.1 <- df[l,1]
        var.2 <- df[l,2]
        z.crit <- (dat_uniq_d[,x.var] == var.1 & dat_uniq_d[,y.var] == var.2)
        z <- mean(na.omit(dat_uniq_d[z.crit,z.var]))
        return(z)
      })
      df$z <- z
      colnames(df) <- c("x","y","z")
      
      if(j == 1) assign(x=paste0("CL2_sp",k,"_K",j),ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(1,17.5)) +
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.x=element_blank(), plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill="none"))
      
      if(j == 2) assign(x=paste0("CL2_sp",k,"_K",j), ggplot(df, aes(x = x, y = (y))) +
                          geom_tile(aes(fill = z)) +
                          scale_fill_viridis(limits = c(1,17.5),breaks=c(1,17)) + 
                          xlab(NULL) +
                          ylab(NULL) +
                          scale_color_viridis() +
                          theme_minimal(base_size=10, base_line_size=0.1) +
                          theme(axis.text.y = element_text(hjust = 20), axis.text.x = element_text(vjust = -1),
                                legend.position="bottom",legend.box="vertical",legend.text = element_text(size=7), 
                                legend.title = element_text(size=0),legend.key.size = unit(0.45, 'cm'),
                                plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
                          annotate("text", x = -0.03, y = 0.36, angle = 90, label = expression(paste( italic(sigma^{P}))), size=4) +
                          coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
                          guides(fill = guide_colorbar(title = " ")))
    }
  }
}

# Productivity gradient - Figure 2 and 3 panel a) and b) 
#################################################################

(6-1)/(2-1) # calculate the slope of the triangle
x1 <- seq(1,2,2/100)
x2 <- seq(2,3,2/100)
y1 <- -4 + (5*x1)
y2 <- 16 + (-5*x2)
datTri <- data.frame(x=c(x1,x2),y=c(y1,y2),co=rep(1,length(c(x1,x2))))

Kgrad1 <- (ggplot(datTri, aes(x, y), show_guide = FALSE) + 
             geom_polygon(alpha=1,colour="darkgray",fill="darkgray") + labs(x = "Efficiency", y = "Mandate") +
             annotate("text", x = 2.1, y = 2.7, col=1, angle = 270, label = str_wrap("Productivity K"), size=3) +
             coord_cartesian(xlim = c(1,4.5), ylim = c(1,6)) +  #, clip = "off"
             theme_void())




#################################################################
### make the figure
#################################################################

# Figure 2 - Correlated case
#################################################################
# save in tiff format
#####################

tiff("Fig2.tiff", units="in", width=6.5, height=6, res=300)
layout <- "
#YY#####ZZ###
#YY#####ZZ###
#YY#####ZZ###
#############
###DDCC#EEFF#
#AADDBB#EEGG#
#AADDBB#EEGG#
#HHIIMJ#NNOOP
#HHIIMJ#NNOOP
#HHIIMJ#NNOOP
#QQRRMK#VVXXP
#QQRRMK#VVXXP
#QQRRMK#VVXXP
"
preyniche1 + predniche1 + predUpniche1 + net1 + net2 + predUpniche2 + predniche2 +
  CL_sp1_K1 + SI_sp1_K1 + plot_spacer() + plot_spacer() + Kgrad1 +  CL_sp2_K1 + SI_sp2_K1 + Kgrad1 +
  CL_sp1_K2 + SI_sp1_K2 + CL_sp2_K2 + SI_sp2_K2 + 
  corr_case + plot_spacer() +
  plot_layout(design = layout)

grid.text("Correlated case",x=unit(0.50, "npc"),y=unit(0.85,"npc"),gp=gpar(fontsize=20)) 
grid.text(expression(paste("Prey niche axis ", italic(nu^{N}) )),x=unit(0.05, "npc"),y=unit(0.895,"npc"),gp=gpar(fontsize=11),rot=90)
grid.text("Prey community",x=unit(0.19, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12, col="#21908CFF"))
grid.text("Predator community",x=unit(0.76, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12, col="#440154FF")) 
grid.text("a)",x=unit(0.07, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12))
grid.text("b)",x=unit(0.62, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12)) 
grid.text("Clusters",x=unit(0.17, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black"))
grid.text("Significance",x=unit(0.32, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Clusters",x=unit(0.71, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Significance",x=unit(0.86, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text(expression(paste(italic(sigma^{N}))),x=unit(0.25, "npc"),y=unit(0.11,"npc"),gp=gpar(fontsize=10, col="black")) 
grid.text(expression(paste(italic(sigma^{N}))),x=unit(0.785, "npc"),y=unit(0.11,"npc"),gp=gpar(fontsize=10, col="black")) 

grid.lines(x=c(95,(1155)), y = c(425,425), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c((1170),(1930))), y = c(425,425), default.units='native' , gp=gpar(col="#440154FF", lty="longdash"))
grid.lines(x=c(95,(1155)), y = (c(1790,1790)), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c((1170),(1930))), y = (c(1790,1790)), default.units='native' , gp=gpar(col="#440154FF", lty="longdash"))

grid.lines(x=(c(95,95)), y = (c(425,1790)), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c(1155,1155)), y = (c(425,1790)), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c(1170,1170)), y = (c(425,1790)), default.units='native' , gp=gpar(col="#440154FF", lty="longdash") )
grid.lines(x=(c(1930,1930)), y = (c(425,1790)), default.units='native' , gp=gpar(col="#440154FF", lty="longdash"))
dev.off()


# Figure 3 - Uncorrelated case
#################################################################

# save in tiff format
#####################

tiff("Fig3.tiff", units="in", width=9, height=6, res=300)
layout <- "
#YY#######ZZ######
#YY#######ZZ######
#YY#######ZZ######
##################
#####DDCC####EEFF#
###AADDBB####EEGG#
###AADDBB####EEGG#
#HHIIJJKKLLM#NNOOP
#HHIIJJKKLLM#NNOOP
#HHIIJJKKLLM#NNOOP
#QQRRSSTTUUM#VVXXP
#QQRRSSTTUUM#VVXXP
#QQRRSSTTUUM#VVXXP
"
prey2niche + pred2niche + predUpniche1 + net1 + net2 + predUpniche2 + predniche2 +  
  CL2_sp1_K1 + SI2_sp1_K1 + plot_spacer() + CL2_sp3_K1 + SI2_sp3_K1 + Kgrad1 +  CL2_sp2_K1 + SI2_sp2_K1 + Kgrad1 +
  CL2_sp1_K2 + SI2_sp1_K2 + plot_spacer() + CL2_sp3_K2 + SI2_sp3_K2 + CL2_sp2_K2 + SI2_sp2_K2 +
  uncorr_case + plot_spacer() +
  plot_layout(design = layout)

grid.text("Uncorrelated case",x=unit(0.50, "npc"),y=unit(0.85,"npc"),gp=gpar(fontsize=20)) 
grid.text(expression(paste("Prey niche axis ", italic(nu^{N}) )),x=unit(0.035, "npc"),y=unit(0.89,"npc"),gp=gpar(fontsize=11),rot=90)
grid.text("Prey community",x=unit(0.13, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12, col="#21908CFF"))
grid.text("Predator community",x=unit(0.79, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12, col="#440154FF")) 
grid.text("a)",x=unit(0.05, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12)) 
grid.text("b)",x=unit(0.69, "npc"),y=unit(0.74,"npc"),gp=gpar(fontsize=12)) 
grid.text("Clusters",x=unit(0.125, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Significance",x=unit(0.235, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Clusters",x=unit(0.455, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Significance",x=unit(0.575, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Clusters",x=unit(0.79, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.text("Significance",x=unit(0.90, "npc"),y=unit(0.03,"npc"),gp=gpar(fontsize=8, col="black")) 
grid.lines(x=c(95,(1800)), y = c(420,420), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c((1815),(2675))), y = c(420,420), default.units='native' , gp=gpar(col="#440154FF", lty="longdash"))
grid.lines(x=c(95,(1800)), y = (c(1787,1787)), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c((1815),(2675))), y = (c(1787,1787)), default.units='native' , gp=gpar(col="#440154FF", lty="longdash"))
grid.text(expression(paste(italic(sigma^{N}))),x=unit(0.18, "npc"),y=unit(0.11,"npc"),gp=gpar(fontsize=10, col="black")) 
grid.text(expression(paste(italic(sigma^{N}))),x=unit(0.52, "npc"),y=unit(0.11,"npc"),gp=gpar(fontsize=10, col="black"))
grid.text(expression(paste(italic(sigma^{N}))),x=unit(0.85, "npc"),y=unit(0.11,"npc"),gp=gpar(fontsize=10, col="black")) 

grid.lines(x=(c(95,95)), y = (c(420,1787)), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c(1800,1800)), y = (c(420,1787)), default.units='native' , gp=gpar(col="#21908CFF", lty="longdash"))
grid.lines(x=(c(1815,1815)), y = (c(420,1787)), default.units='native' , gp=gpar(col="#440154FF", lty="longdash") )
grid.lines(x=(c(2675,2675)), y = (c(420,1787)), default.units='native' , gp=gpar(col="#440154FF", lty="longdash"))
dev.off()
