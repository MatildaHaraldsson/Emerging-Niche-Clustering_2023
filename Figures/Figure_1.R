#######################################################################################
# Script for making Figure 1 in Haraldsson and Thébault 2023, Ecology Letters
#######################################################################################

rm(list=ls())

# libraries
library(ggplot2)
library(viridis)
library(tidyverse)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)
library(RGraphics)
library(CircStats)
library(gtable)

source("Library_Emerging_niche_clustering.R")

# load necessary data
(load("datClust_spnr2_t500_852226.RData")) # simulation and cluster results from simulation 852226


#################################################################
### a) 
#################################################################

# Competition and predation kernels 
#################################################################

mu_N <- seq(0,1,1/100) # for generating alpha values
a_N <- alpha_N(mu=mu_N,sigma_sp=0.10,nr_sp=length(mu_N),Lmax=1,P=4,Sigma_term=1) # generate alpha values for plotting competition kernel
a_P <- beta_P(mu=mu_N,mu_pred=mu_N,sigma_sp=0.20,sigma_P=0.20,nr_sp=length(mu_N),nr_P=length(mu_N),Lmax=1,P=4,Sigma_term=1) # generate beta values for plotting predation kernel

colnames(a_N) <- round(mu_N,digits=3)
rownames(a_N) <- round(mu_N,digits=3)
colnames(a_P) <- round(mu_N,digits=3)
rownames(a_P) <- round(mu_N,digits=3)

dat <- as.data.frame(cbind(a=a_N[,"0.4"],mu=mu_N))  # target species at mu = 0.4 in the competition matrix
datP <- as.data.frame(cbind(a=a_P[,"0.4"],mu=mu_N)) # target species at mu = 0.4 in the predation matrix
datPointN <- as.data.frame(cbind(a=rep(0,length(sort(runif(20,0,1)))),mu=seq(0,1,length.out=20)))

# prey niche axis
preyniche <- (ggplot(NULL,aes(y=a,x=mu)) +
  geom_line(data=dat, colour="darkgray") +
  geom_line(data=datP, colour="darkgray", linetype="dashed") +
  geom_point(data=datPointN, shape=21, size=1, color="#21908CFF") + 
  geom_point(aes(x=0.4,y=0), size=3, color="#21908CFF") + 
  coord_cartesian(xlim = c(0,1),ylim = c(0,(1.2*1.2)), clip = "off") +
  scale_x_continuous(breaks=seq(0,1,0.5),limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab(expression(paste("Prey niche axis (", italic(mu^{N}),")  "))) +
  ylab(NULL) +
  geom_segment(aes(x=0.4,y=0.6,xend=0.485,yend=0.6), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.485,y=0.6,xend=0.4,yend=0.6), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.4,y=-0.4,xend=0.4,yend=0),arrow=arrow(angle=9,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0.4,y=0.4,xend=0.595,yend=0.4), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_segment(aes(x=0.595,y=0.4,xend=0.4,yend=0.4), arrow=arrow(angle=20,length=unit(0.05,"inches"))) +
  geom_line(data=data.frame(cbind(mu=c(0.4,0.4,0.4),a=c(0,0.6,1.2))), linetype="dotdash",color="darkgray") +
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
)

# Simplified network
#################################################################
datA <- data.frame(x=c(0,3,6),y=c(2.95,2.95,2.95))
datB <- data.frame(x=c(0,3,6),y=c(1,1,1))
datC <- data.frame(x=c(0,3,6),y=c(-0.95,-0.95,-0.95))

net <- ggplot(NULL,aes(x=x, y=y)) + 
  theme_void() +
  scale_fill_manual(values=c("darkgray","#21908CFF","darkgray")) +
  
  geom_segment(aes(x=0,y=1,xend=0,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=1,xend=3,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=1,xend=6,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=3,y=1,xend=0,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=1,xend=3,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=1,xend=6,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=6,y=1,xend=0,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=1,xend=3,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=1,xend=6,yend=2.33),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=0,y=-0.95,xend=0,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=-0.95,xend=3,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=0,y=-0.95,xend=6,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=3,y=-0.95,xend=0,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=-0.95,xend=3,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=3,y=-0.95,xend=6,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_segment(aes(x=6,y=-0.95,xend=0,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=-0.95,xend=3,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  geom_segment(aes(x=6,y=-0.95,xend=6,yend=0.38),arrow=arrow(angle=7,length=unit(0.10,"inches"),type="closed")) +
  
  geom_point(data=datC,shape=21,aes(fill="datC"),size=8,show.legend = FALSE) +
  geom_point(data=datA,shape=21,aes(fill="datA"),size=8,show.legend = FALSE) +
  geom_point(data=datB,shape=21,aes(fill="datB"),size=8,show.legend = FALSE) +
  
  coord_cartesian(xlim = c(-0.80,5.2), ylim = c(-0,3), clip = "off")

#################################################################
### b)
#################################################################

# Identified clusters from one simulation
#################################################################

#  Parameters for simulation nr 852226:
#  K <- 5
#  c <- 0.6
#  sN <- 0.30
#  sP <- 0.30
#  r <- 0.5926975
#  Q <- 4
#  n <- 200

  data_plot <- data.frame(datNoZeroClust)
  data_plot <- (cbind(data_plot,N.freq=(data_plot[,"N.abu"]/sum(data_plot[,"N.abu"])),N.percent=(data_plot[,"N.abu"]/sum(data_plot[,"N.abu"])*100)))  

  niche_plot <- ggplot(data_plot, aes(N.mu,N.abu)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, (max(data_plot$N.abu)*1.1))) +
    scale_y_continuous(name="Abundance", expand = c(0, 0), breaks=c(0, .4, .8)) +
    scale_x_continuous(expand = c(0, 0),breaks=c(0,0.50,1)) +
    geom_point(size=0.5,color=1,alpha=0) +
    xlab(NULL) + # xlab(expression(paste( italic(?^{N})))) +
    geom_linerange(x=data_plot$N.mu,ymin=rep(0,dim(data_plot)[1]),ymax=data_plot$N.abu) +
    theme_minimal(base_size=10, base_line_size=0.1)

# Eigenvalues for competition case
#################################################################
  compet<-function(muN,sN,Q){
    a<-NULL
    for (i in 1:length(muN)){
      for (j in 1:length(muN))
        a<-c(a,exp(-(min(c(abs(muN[i]-muN[j]),(1-abs(muN[i]-muN[j]))))/(sN))^Q))
    }
    return(matrix(a,nrow=length(muN)))
  }
  
  matunit<-function(n){
    u<-matrix(1,nrow=n,ncol=n)
    for (i in 2:n){
      for (j in 2:n)
        u[j,i]<-complex(real=cos(2*pi*(i-1)*(j-1)/n),imaginary=sin(2*pi*(i-1)*(j-1)/n))
    }
    return(u)
  }
  
  K <- 5
  sN <- 0.30
  r <- 0.5926975
  Q <- 4
  n <- 200
  sN <- seq(0.1,0.3,by=0.01)
  sN <- 0.30
  valmax <- NULL
  nbclust <- NULL
  n <- 200

  for (i in 1:length(sN)) {
    muN<- seq(from=0,by=(1/n),length.out=n) # runif(n,0,1) #  # higher resolution in image
    a<-compet(muN,sN[i],Q)
    Nstar<-K/sum(a[1,])
    U<-matunit(n)
    eigenvaltot<-NULL
    Z<-a
    for (j in 1:(n/2+1)){
      eigenvaltot<-c(eigenvaltot,-1*(Z[1,]%*%U[,j]))
    }
    index<-which(Re(eigenvaltot)==max(Re(eigenvaltot)))
    valmax<-c(valmax,Re(eigenvaltot[index]))
    nbclust<-c(nbclust,index-1)
  }
  dat_eig <- data.frame(cbind(muN,eigen=Re(U[,nbclust[length(sN)]+1])))
  
  eigen_plot <- ggplot(dat_eig, aes(muN,eigen)) +
    coord_cartesian(xlim = c(0, 1), ylim = c((min(dat_eig$eigen)), (max(dat_eig$eigen)*1.1))) +
    scale_y_continuous(name="Eigenvector",expand = c(0, 0),breaks=c(-1, 0, 1)) +
    scale_x_continuous(expand = c(0, 0),breaks=c(0,0.50,1)) +
    geom_point(size=0.5,color=1,alpha=0.5) +
    xlab(expression(paste( italic(mu^{N})))) +
    theme_minimal(base_size=10, base_line_size=0.1)

#################################################################
### c)
#################################################################

# Relative importance of predation ###
#################################################################
# The calculations for the analytical results (below) are time consuming.
#     Alternative 1) load the output directly
#     Alternative 2) calculations analytical results via function calcAnalytical() with RelImpPred = 0.01, 1 and 100
  
#     Alternative 1) load the results
  
(load("nbclustmat_Q4_aR0_01.RData"))
eigenDat01 <- as.data.frame(cbind(nbclust=c(nbclustmat),
                                  sP=rep(as.numeric(rownames(nbclustmat)),ncol(nbclustmat)),
                                  sN=rep(as.numeric(colnames(nbclustmat)),each=nrow(nbclustmat))))
(load("nbclustmat_Q4_aR1.RData"))
eigenDat1 <- as.data.frame(cbind(nbclust=c(nbclustmat),
                                 sP=rep(as.numeric(rownames(nbclustmat)),ncol(nbclustmat)),
                                 sN=rep(as.numeric(colnames(nbclustmat)),each=nrow(nbclustmat))))
(load("nbclustmat_Q4_aR100.RData"))
eigenDat100 <- as.data.frame(cbind(nbclust=c(nbclustmat),
                                   sP=rep(as.numeric(rownames(nbclustmat)),ncol(nbclustmat)),
                                   sN=rep(as.numeric(colnames(nbclustmat)),each=nrow(nbclustmat))))
eigenDat <- (cbind(rbind(eigenDat01,eigenDat1,eigenDat100),relPred=c(rep(0.01,676),rep(1,676),rep(100,676))))


#     Alternative 2) calculate the results

calcAnalytical <- function(RelImpPred) {
# Model Chesson & Kuang

# Expected number of clusters as function of sigma N and sigma P for different levels of Importance of Predation 
n <- 200
Q <- 4      
rP <- 1
rR <- 1
alphaP <- 1  
alphaR <- RelImpPred       # For Figure 1, the relative importance of predation = [0.01, 1, and 100].
sN <- seq(0.05,0.3,by=0.01)
sP <- seq(0.05,0.3,by=0.01)

nbclustmat<-matrix(0,nrow=length(sP),ncol=length(sN))
rownames(nbclustmat)<-sP
colnames(nbclustmat)<-sN
muN<-seq(from=0,by=(1/n),length.out=n)
U<-matunit(n)

for(nN in 1:length(sN)) {
  for(nP in 1:length(sP)) {
    theta<-compet(muN,sN[nN],Q)
    gamma<-compet(muN,sP[nP],Q)
    valmax<-NULL
    nbclust<-NULL
    for (i in 1:length(alphaR)){
      eigenvaltot<-NULL
      beta<-1/(rR*alphaR[i])*theta + 1/(rP*alphaP)*gamma
      for (j in 1:(n/2+1)){
        eigenvaltot<-c(eigenvaltot,-1*(beta[1,]%*%U[,j]))
      }
      index<-which(Re(eigenvaltot)==max(Re(eigenvaltot)))
      valmax<-c(valmax,Re(eigenvaltot[index]))
      nbclust<-c(nbclust,index-1)
    }
    nbclustmat[nP,nN] <- nbclust
  }
}

return(nbclustmat)

} # end of function

# Unmark if calculating the analytical results
# anaRes001 <- calcAnalytical(RelImpPred=0.01) 
# anaRes1 <- calcAnalytical(RelImpPred=1)  
# anaRes100 <- calcAnalytical(RelImpPred=100)  

marg1 <- marg3 <- 0.42
marg2 <- marg4 <- 0.00

eigenMy <- ggplot(eigenDat, aes(x = sN, y = (sP))) +
  geom_tile(aes(fill = nbclust)) +
  scale_fill_viridis(limits = c(2,16.2),breaks=c(2,16)) + 
  xlab(expression(paste( italic(sigma^{N})))) +
  ylab(expression(paste( italic(sigma^{P})))) +
  scale_color_viridis() +
  theme_minimal(base_size=10, base_line_size=0.1) +
  theme(axis.text.y = element_text(hjust = 20), axis.text.x = element_text(vjust = -1),
        legend.position="bottom",legend.box="vertical",legend.text = element_text(size=7),
        legend.title = element_text(size=0),legend.key.size = unit(0.45, 'cm'),
        plot.margin = unit(c(marg1,marg2,marg3,marg4), "cm")) +
  coord_cartesian(xlim = c(0.05,0.30), ylim = c(0.05,0.30), clip = "off") +
  guides(fill = guide_colorbar(title = " ")) +
  facet_wrap(~ relPred, ncol = 1) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Productivity gradient
#################################################################
(6-1)/(2-1) # calculate the slope of the triangle
x1 <- seq(1,2,2/100)
x2 <- seq(2,3,2/100)
y1 <- -4 + (5*x1)
y2 <- 16 + (-5*x2)
datTri <- data.frame(x=c(x1,x2),y=c(y1,y2),co=rep(1,length(c(x1,x2))))

Kgrad <- (ggplot(datTri, aes(x, y), show_guide = FALSE) + 
            geom_polygon(alpha=1,colour="darkgray",fill="darkgray") + labs(x = "Efficiency", y = "Mandate") +
          annotate("text", x = 3.8, y = 3, col=1, angle = 270, label = str_wrap("Relative importance of predation"), size=3) +
        coord_cartesian(xlim = c(-1,4.5), ylim = c(1,6)) +
        theme_void())

#################################################################
### make the figure
#################################################################

# save in tiff format
#####################
tiff("Fig1.tiff", units="in", width=5.016, height=4.37, res=300)

layout <- "
#AABBEEH
#AABBEEH
#AABBEEH
#CCCCEEH
#CCCCEEH
#CCCCEEH
#DDDDEEH
#DDDDEEH
#DDDDEEH
"
preyniche + net + 
  niche_plot + eigen_plot + 
  eigenMy + Kgrad +
  plot_layout(design = layout)

grid.text("a)",x=unit(0.07, "npc"),y=unit(0.98,"npc"))
grid.text("b)",x=unit(0.07, "npc"),y=unit(0.70,"npc")) 
grid.text("c)",x=unit(0.66, "npc"),y=unit(0.98,"npc")) 
grid.text(expression(italic(sigma^{N})),x=unit(0.31, "npc"),y=unit(0.89,"npc")) 
grid.text(expression(italic(sigma^{P})),x=unit(0.33, "npc"),y=unit(0.86,"npc")) 
grid.text("Clusters",x=unit(0.79, "npc"),y=unit(0.05,"npc"),gp=gpar(fontsize=7.5,col="black")) 

dev.off()

# save in pdf format
#####################

dir_plot <- paste0(getwd(),"/Fig1.pdf") # save the plot to the current work directory
pdf(file = dir_plot,   # filename and chosen directory
    width = 4.5,       # the width of the plot in inches
    height = 4.37      # the height of the plot in inches
)

layout <- "
AABBEEH
AABBEEH
AABBEEH
CCCCEEH
CCCCEEH
CCCCEEH
DDDDEEH
DDDDEEH
DDDDEEH
"
preyniche + net + 
  niche_plot + eigen_plot + 
  eigenMy + Kgrad +
  plot_layout(design = layout)

grid.text("a)",x=unit(0.02, "npc"),y=unit(0.98,"npc")) 
grid.text("b)",x=unit(0.02, "npc"),y=unit(0.70,"npc")) 
grid.text("c)",x=unit(0.61, "npc"),y=unit(0.98,"npc")) 
grid.text(expression(italic(sigma^{N})),x=unit(0.26, "npc"),y=unit(0.89,"npc"))
grid.text(expression(italic(sigma^{P})),x=unit(0.28, "npc"),y=unit(0.86,"npc")) 
grid.text("Clusters",x=unit(0.77, "npc"),y=unit(0.05,"npc"),gp=gpar(fontsize=7.5,col="black")) 

dev.off()