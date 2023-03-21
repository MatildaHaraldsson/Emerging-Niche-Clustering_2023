# Library with functions for simulating predator prey systems with a single or double niches axes for prey
# as in Haraldsson & Thébault 2023
#######################################################################################

# -----------------------------------------------------------------------

# demanded packages
require("deSolve")   # for ode simulation
require("ggplot2")   # for nicheplot
require("patchwork") # for nicheplot

# -----------------------------------------------------------------------
# Model_Prey_Pred() - ode model: competing prey and predators
# -----------------------------------------------------------------------

Model_Prey_Pred <- function(t,y,parms){ 

  nr_sp <- parms[1]               # nr competing species
  nr_P  <- parms[2]               # nr predators
  r_par <- parms[3]               # maximum per capita growth rate of species i
  K_par <- parms[4:(3+nr_sp)]     # carrying capacity of species i, assumed to be the same for all species
  m_par_pred <- parms[(4+nr_sp)]  # mortality of predators
  e_par_pred <- parms[(5+nr_sp)]  # conversion factor from one trophic level to the next
  c_par <- parms[(6+nr_sp)]       # maximum attack rate
  
  a_par <- parms[(7+nr_sp):(6+nr_sp+nr_sp^2)]                            # competition coefficients alpha (community matrix)
  a_par_pred <- parms[(7+nr_sp+nr_sp^2):(6+nr_sp+nr_sp^2+nr_P*nr_sp)]    # attack rate beta
  
  N <- y[1:nr_sp]                 # denstiy of competing species
  P <- y[(nr_sp+1):(nr_sp+nr_P)]  # density of predators  
  
  yprime<-0*y
  
  N[N < 0] <- 0
  P[P < 0] <- 0
  
  A <- matrix(a_par,nrow=nr_sp,ncol=nr_sp,byrow=TRUE)
  A_pred <- matrix(a_par_pred,nrow=nr_sp,ncol=nr_P,byrow=TRUE)
  
  dN <- r_par * N * (K_par - (A %*% N)) / K_par - c_par * (A_pred %*% P * N)
  dP <- e_par_pred * c_par * (t(A_pred) %*% N * P) - (P * m_par_pred)
  
  yprime[1:nr_sp] <- dN
  yprime[(nr_sp+1):(nr_sp+nr_P)] <- dP
  
  return(list(yprime,c()))
}


# -----------------------------------------------------------------------
# Model_Comp - ode model: competing species without predation
# -----------------------------------------------------------------------

Model_Comp <- function(t,y,parms){ 
  nr_sp <- parms[1]               # nr competing species
  r_par <- parms[3]               # maximum per capita growth rate of species i
  K_par <- parms[4:(3+nr_sp)]     # carrying capacity of species i # assumed to be the same for all species in Scheffer & van Nes 2006
  
  a_par <- parms[(7+nr_sp):(6+nr_sp+nr_sp^2)]                            # competition coefficients alpha (community matrix)
  
  N <- y[1:nr_sp]                 # denstiy of competing species
  
  yprime<-0*y[1:nr_sp]
  
  N[N < 0] <- 0
  
  A <- matrix(a_par,nrow=nr_sp,ncol=nr_sp,byrow=TRUE)
  
  dN <- r_par * N * (K_par - (A %*% N)) / K_par
  
  yprime[1:nr_sp] <- dN
  
  return(list(yprime,c()))
}

# -----------------------------------------------------------------------
# dwrp_flex() - estimates the wrapped density function
# -----------------------------------------------------------------------

dwrp_flex <- function (theta, mu, rho, P = 2, sigma_term = 1, sd = 1, acc = 1e-05, tol = acc) {
  # Modified from the dwrpnorm() function in the CircStats library
  # Including a flexible form of the density function depending on parameter P, 
  # following Pigolotti et al. 2010,  Theor. Ecol. 3: 89-96. 
  
  # P determines the shape of the distribution. When P < 2 the kernels are more peaked around y ~ 0
  # and when P > 2 they become more box-like
  #     P = 1 exponential distribution
  #     P = 2 normal distribution
  #     P = 3 or 4 box distribution
  # sigma_term : adjusts the shape, see Pigolotti et al. 2010 (page 91...)
  
  if (missing(rho)) {
    rho <- exp(-sd^2/2)
  }
  if (rho < 0 | rho > 1) 
    stop("rho must be between 0 and 1")
  var <- -2 * log(rho)
  term <- function(theta, mu, var, k, P, sigma_term) {
    exp(- abs(( (theta - mu + 2 * pi * k)^P)/(sigma_term * sqrt(var)^P)))
  }
  k <- 0
  Next <- term(theta, mu, var, k, P, sigma_term)
  delta <- 1
  while (delta > tol) {
    k <- k + 1
    Last <- Next
    Next <- Last + term(theta, mu, var, k, P, sigma_term) + term(theta, mu, var, -k, P, sigma_term)
    delta <- abs(Next - Last)
  }
  Next
}

# -----------------------------------------------------------------------
# alpha_N() - computes competition coefficient (alpha) scaling the effect of competitor j on i 
# -----------------------------------------------------------------------

alpha_N <- function(mu,sigma_sp,nr_sp,Lmax=1,P=4,Sigma_term=1) {
  # assumes the width of the niche follows the distribution defined in the dwrp_flex() function
  
  # mu         = trait values along niche axis
  # sigma_sp   = niche width
  # nr_sp      = number of competing species
  # Lmax       = length of niche axis
  # P          = parameter determining the shape of the stretched exponential function
  # Sigma_term = 1
  
  # returns the community matrix A with competition coefficients
  
  if(NROW(sigma_sp) == 1) sigma_sp <- rep(sigma_sp,nr_sp)
#  source("dwrpnorm_flex.R")
  
  A <- matrix(NA,nr_sp,nr_sp)
  for (i in 1:nr_sp) {   # competing species i
    for (j in 1:nr_sp) { # competing species j
      muimoinsj_rad <- (mu[i]-mu[j])*2*pi/Lmax   # trait difference between competing species i and j expressed in radians
      sigma_compi <- sigma_sp[i]*2*pi/Lmax       # niche width of species i expressed in radians
      a <- dwrp_flex(0,muimoinsj_rad,sd=sigma_compi,P=P,sigma_term=Sigma_term)
      A[i,j] <- a
    }
  }
  return(A)
}

# -----------------------------------------------------------------------
# beta_P() - computes the attack rate (beta) of consumer j on their prey i
# -----------------------------------------------------------------------

beta_P <- function(mu,mu_pred,sigma_sp,sigma_P,nr_sp,nr_P,Lmax=1,P=4,Sigma_term=1) {
  # assumes the width of the niche follows the distribution defined in the dwrp_flex() function
  
  # mu / mu_pred         = trait values along niche axis for prey or predators
  # sigma_sp / sigma_P   = niche width
  # nr_sp / nr_P         = number of competing species
  # Lmax                 = length of niche axis
  # P                    = parameter determining the shape of the stretched exponential function
  # Sigma_term           = 1
  
  # returns the matrix B with attack rates
  
  if(NROW(sigma_sp) == 1) sigma_sp <- rep(sigma_sp,nr_sp)
  if(NROW(sigma_P) == 1) sigma_P <- rep(sigma_P,nr_P)
#  source("dwrpnorm_flex.R")
  
  B <- matrix(NA,nr_sp,nr_P)
  for (i in 1:nr_sp) {    # competing prey i
    for (j in 1:nr_P) {   # predator j
      muimoinsj_rad <- (mu[i]-mu_pred[j])*2*pi/Lmax   # trait difference between prey i and predator j expressed in radians
      sigma_predj <- sigma_P[j]*2*pi/Lmax  
      a <- dwrp_flex(0,muimoinsj_rad,sd=sigma_predj,P=P,sigma_term=Sigma_term)
      B[i,j] <- a
    } 
  } 
  return(B)
}

# -----------------------------------------------------------------------
# MadedatForClust() - takes ode simulation data and format it for cluster analysis using KmeansGap() 
# -----------------------------------------------------------------------

MakedatForClust <- function(res,muList,niche=NULL,treshold=NULL,TExclude=TRUE) {
  # Input:
  #       res       = output from ode simulation
  #       muList    = list with trait values (mu_N, mu_P, and mu_N)
  #       treshold  = extinction treshold        parms - parameter settings used in lsoda
  #       TExclude  = TRUE/FALSE, excludes extinct species from dataset if TRUE
  # Output:
  #       dat       = list with data formatted for analysis with KmeansGap()
  
  if(NCOL(res) == (nr_spN+nr_spP+1)) {
    if(niche == "corr") spp <- c("N","P")
    if(niche == "uncorr") spp <- c("N","P","M")
  } else {
    spp <- c("N")
  }
  dat <- sapply(spp, function(sp) {
    if(sp == "N")  trait <- muList[["mu_N"]]  
    if(sp == "P")  trait <- muList[["mu_P"]]
    if(sp == "M") { sp <- "N" ; trait <- muList[["nu_N"]] }
    N <- res[dim(res)[1],grep(sp,colnames(res))]
    dat <- cbind(N,trait)
    if(!is.null(treshold)) dat[which(dat[,"N"] < treshold),"N"] <- 0 # applies extinction threshold and equal to zero
    if(TExclude) dat <- dat[which(dat[,"N"] > 0),] # if TRUE, exclude all zeros from dataset
    return(dat)
  }, simplify = FALSE,USE.NAMES = TRUE)
  return(dat)
}


# -----------------------------------------------------------------------
# runKmeansGap() - runs a modified version of KmeansGap, 
# origingally by D'Andrea et al (2019) Generalizing clusters of similar species as a signature of coexistence under competition, PLoS Comp. Biol. 15(1): e1006688.
# https://github.com/rafaeldandrea/Clustering-metric
# the function includes the arguments used in the simulations in Haraldsson & Thébault 2023
# and returns the results used for further analysis and figures
# -----------------------------------------------------------------------

runKmeansGap <- function(datForClust) {
  # datClust        = dataformat for KmeansGap, abundances "N" and trait values "trait". Can be generated though datForClust() function.
  #
  # Output: clust, a list for each community that was entered via the data "datForClust" (N, P and/or M) containing two lists (clustSum, clusterID) of main results from the KmeansGap().
  # list - clustSum:
  #     khat        = number of cluster 
  #     maxGap      = Gap statistic, maximum value of the gap index across all number of clusters tested in the observed community 
  #     Zscore      = Defined as (maxgap-mean(maxnullgap))/sd(maxnullgap), see D'Andrea et al. 2019.
  #     Pvalue      = P-value
  # list - clusterID:
  #                 = data.frame containing original abundances (abu), trait value (mu) and cluster identity (clustID) for each specimen. Contains only speciment with abundances > 0.
  
  spp <- names(datForClust)
  
  clustRes <- sapply(spp, function(sp) {
    
    Dat <- data.frame(datForClust[[sp]])
    
    if(NROW(Dat) > 1) {
      clust <- KmeansGapMod(  
        dat=Dat,
        nozeros=TRUE,
        multiD=FALSE,
        nullmodel="shuffle",
        numnulls=100,
        mink=1,
        maxk=20,
        nstartingpoints=100,
        weighting="yan",
        shortcut=FALSE,
        internal.peaks=FALSE,
        plot=FALSE,
        bands=TRUE,
        plotquant90=TRUE,
        verbose=TRUE)
      
      cID <- clust$clusterID[,c(clust$khat+2)]
      DatClust <- merge(Dat,cID,by="row.names")
      rownames(DatClust) <- DatClust[,"Row.names"]
      DatClust <- DatClust[,c(2:4)]
      colnames(DatClust) <- c("abu","mu","clustID")
      DatClust <- DatClust[(order(DatClust$mu)),]
    }
    clustSum <- data.frame(khat=if(NROW(Dat) > 1) clust$khat else NA,
                           maxgap=if(NROW(Dat) > 1) clust$maxgap else NA,
                           zscore=if(NROW(Dat) > 1) clust$z.score else NA,
                           pvalue=if(NROW(Dat) > 1) clust$p.value else NA)
    
    listRes <- list(clustSum=clustSum,
                    clusterID=if(NROW(Dat) > 1) DatClust else NA)
    return(listRes)
  }, simplify = FALSE,USE.NAMES = TRUE)
  
  return(clustRes) 
}

# -----------------------------------------------------------------------
# nicheplotClust() - plots the abundance of species along its niche axis with the clusters indicated in alternating colors
# -----------------------------------------------------------------------

nicheplotClust <- function(datClust,TL=2) {
  # Uses the output from the runKmeansGap() function
  # arrange data to plot the abundance of a trait along the niche axis 
  #
  # datClust        = the output from the runKmeansGap() function as a list for each species. Contains:
  #       abu       = abundances of N or P at a given timestep (typically last time step in a lsoda() simulation)
  #       mu        = trait value
  #       clusterID = cluster identificatio number
  
  spp <- names(datClust)
  dat <- datClust
  # set parameters for plotting depending on dataset
  if(length(spp) == 1) {
    if(spp[1] == "N") { 
      nr_sp <- nr_sp2 <- dim(dat[[spp[1]]]$clusterID)[1] 
      p_title <- c("Prey",NA,NA)
      x_title <- c(expression(paste("Prey niche axis ", italic(mu^{N}))),NA,NA)
      colPlot <- c("#21908CFF")
      if(TL == 1) { p_title <- c("Competing community",NA,NA)
      x_title <- c(expression(paste("Niche axis ", italic(mu^{N}))),NA,NA)
      colPlot <- c("#FDE725FF")
      }
    }
    if(spp[1] == "P") {
      nr_p  <- nr_sp2 <- dim(dat[[spp[1]]]$clusterID)[1] 
      p_title <- c("Predators",NA)
      x_title <- c(expression(paste("Predator niche axis ", italic(mu^{P}))),NA,NA)
      colPlot <- c("#453781FF")
    }
    layout = "
    A
    "} 
  if(length(spp) == 2) {
    if(spp[1] == "N" & spp[2] == "P") {
      nr_sp <- dim(dat[[spp[1]]]$clusterID)[1] 
      nr_p <- dim(dat[[spp[2]]]$clusterID)[1]
      nr_sp2 <- c(nr_sp,nr_p)
      p_title <- c("Prey","Predator")
      x_title <- c(expression(paste("Prey niche axis ", italic(mu^{N}))), expression(paste("Predator niche axis ", italic(mu^{P}))),NA)
      colPlot <- c("#21908CFF","#453781FF")
    }
    if(spp[1] == "N" & spp[2] == "M"){
      nr_sp <- dim(dat[[spp[1]]]$clusterID)[1] 
      nr_p <- dim(dat[[spp[2]]]$clusterID)[1]
      nr_sp2 <- c(nr_sp,nr_p)
      p_title <- c("Prey on competition niche (mu)","Prey on predator niche (nu)")
      x_title <- c(expression(paste("Prey niche axis ", italic(mu^{N}))),expression(paste("Prey niche axis ", italic(nu^{N}))),NA)
      colPlot <- c("#21908CFF","#73D055FF")
    }
    layout = "
    A
    B
    "}
  
  if(length(spp) == 3) {
    nr_sp <- dim(dat[[spp[1]]]$clusterID)[1] 
    nr_p <- dim(dat[[spp[2]]]$clusterID)[1]
    nr_m <- dim(dat[[spp[3]]]$clusterID)[1]
    nr_sp2 <- c(nr_sp,nr_p,nr_m)
    p_title <- c("Prey on competition niche (mu)","Predator","Prey on predator niche (nu)")
    x_title <- c(expression(paste("Prey niche axis ", italic(mu^{N}))),expression(paste("Predator niche axis ", italic(mu^{P}))),expression(paste("Prey niche axis ", italic(nu^{N}))))
    colPlot <- c("#21908CFF","#453781FF","#73D055FF")
    layout = "
    A
    B
    C
    "}
  plotit <- lapply(c(1:length(spp)), function(p) {
    color = dat[[p]]$clusterID[,3]  
    dat[[p]]$clusterID <- dat[[p]]$clusterID[order(dat[[p]]$clusterID$mu),]
    nrClust <- unique(dat[[p]]$clusterID[,3])
    colSeq <- rep(c(colPlot[p],"#000000"),ceiling(max(nrClust)/2))[1:max(nrClust)]
    (color <- as.data.frame(color))
    colrange <- unique(color)
    for(l in 1:max(colrange)) color <- replace(color,color == colrange[l,],colSeq[l])
    pchRange <- replace(color,color == colPlot[p],1)
    pchRange <- replace(pchRange,color == "#000000",0.2)
    colnames(pchRange) <- "pch"
    x <- dat[[p]]$clusterID[,"mu"]
    y <- (dat[[p]]$clusterID[,"abu"])
    plotit <- data.frame(x,y,color,pchRange)  
    return(plotit)
  })
  for(np in 1:length(spp)) {
    data_plot <- data.frame(plotit[[np]])
    assign(paste0("p",np),value=ggplot(data_plot, aes(x,y)) +
             ggtitle(p_title[np]) +
             coord_cartesian(xlim = c(0, 1), ylim = c(0,(max(data_plot$y)*1.01))) +
             scale_y_continuous(name="Abundance", expand = c(0, 0), breaks=c(0,((signif(max(data_plot$y),digits=2))/2),(signif(max(data_plot$y),digits=2)) )) +
             scale_x_continuous(expand = c(0, 0),name=x_title[np],breaks=c(0,0.50,1)) +
             geom_point(size=1.5,color=data_plot$color,alpha=as.numeric(data_plot$pch)) +
             geom_linerange(x=data_plot$x,ymin=rep(0,dim(data_plot)[1]),ymax=data_plot$y,color=data_plot$color) +
             theme_classic(base_size=10))
  }
  windows(width=3.5, height=(3*length(spp))) # plot
  if(length(spp) == 1) print(p1 + plot_layout(design = layout))
  if(length(spp) == 2) print(p1 + p2 + plot_layout(design = layout))
  if(length(spp) == 3) print(p1 + p2 + p3 + plot_layout(design = layout))
}

# -----------------------------------------------------------------------