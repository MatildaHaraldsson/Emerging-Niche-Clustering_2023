#######################################################################################
# Main script for ode-simulations, cluster analysis and richness calculations
# Haraldsson and Thébault 2023
#
# Simulation and analysis for:
# a system with prey and predators with uncorrelated niche axes responsible for competition between prey and predator interaction
# a system with prey and predators with correlated niche axes responsible for competition between prey and predator interaction
# a system with competing species only 
#######################################################################################

library("deSolve")                           # for ode simulation
library("ggplot2")                           # for nicheplot() function
library("patchwork")                         # for Nicheplot() function
library("plyr")                              # for KmeansGap() function
library("pracma")                            # for KmeansGap() function

source("Library_Emerging_niche_clustering.R") # library with functions
source("KmeansGapMod.R")                  # modified version of KmeansGap() function
(load("parameters.R"))                       # parameters for the ode simulations

# -----------------------------------------------------------------------

# NOTE: set simulation ID, between 1 and 874.000
i <- 87353

# -----------------------------------------------------------------------
# parameters
# -----------------------------------------------------------------------

# fixed parameters
times <- seq(from=0,by=0.01,to=500)
nr_spN <- 200                               # number of competing species
nr_spP <- 200                               # number of predator species
lmax <- 1                                   # niche length

# parameters with added variation within a given range
r <- parameters[i,"r"]
N0 <- parameters[i,"N0"]
m_pred <- parameters[i,"m_pred"]
e_pred <- parameters[i,"e_pred"]
P0 <- parameters[i,"P0"]

# parameters specific for simulation
K <- parameters[i,"K"]                      # carrying capacity
K <- rep(K,nr_spN)
sigma_spN <- parameters[i,"sigma_spN"]      # niche width of competing prey 
sigma_spP <- parameters[i,"sigma_spP"]      # niche width of predators
c_ <- parameters[i,"c_"]                    # maximum attack rate
para_sim <- c(nSim=i,K=K[1],c_=c_,sigma_spN=sigma_spN,sigma_spP=sigma_spP,r=r,m_pred=m_pred,e_pred=e_pred,N0=N0,P0=P0) # save vector of parameters for later analysis
para_sim_comp <- c(nSim=i,K=K[1],c_=c_,sigma_spN=sigma_spN,r=r,N0=N0) # save vector of parameters for later analysis


# random parameters - random trait values along a correlated (mu) or uncorrelated (mu and nu) niche axes
set.seed(parameters[i,"mu_N"])
mu_N <- runif(nr_spN,0,lmax)                # competing species/prey niche axis (mu N)
set.seed(parameters[i,"mu_P"])
mu_P <- runif(nr_spP,0,lmax)                # predator niche axis (mu P)
set.seed(parameters[i,"nu_N"])
nu_N <- runif(nr_spN,0,lmax)                # competing prey niche axis (nu N) (for uncorrelated case) 
muList <- list(mu_N=mu_N,mu_P=mu_P,nu_N=nu_N)

# correlated case
a_N <- as.vector(t(alpha_N(mu=mu_N,sigma_sp=sigma_spN,nr_sp=nr_spN,Lmax=1,P=4,Sigma_term=1)))
a_P_corr <- as.vector(t(beta_P(mu=mu_N,mu_pred=mu_P,sigma_sp=sigma_spN,sigma_P=sigma_spP,nr_sp=nr_spN,nr_P=nr_spP,Lmax=lmax,P=4,Sigma_term=1)))
para_corr <- c(nr_sp=nr_spN,nr_P=nr_spP,r=r,K=K,m_pred=m_pred,e_pred=e_pred,c_=c_,a_N=a_N,a_P=a_P_corr) # save vector of parameters for ode simulation

# uncorrelated case
a_P_uncorr <- as.vector(t(beta_P(mu=nu_N,mu_pred=mu_P,sigma_sp=sigma_spN,sigma_P=sigma_spP,nr_sp=nr_spN,nr_P=nr_spP,Lmax=lmax,P=4,Sigma_term=1)))
para_uncorr <- c(nr_sp=nr_spN,nr_P=nr_spP,r=r,K=K,m_pred=m_pred,e_pred=e_pred,c_=c_,a_N=a_N,a_P=a_P_uncorr) # save vector of parameters for ode simulation

# -----------------------------------------------------------------------
# run ode-simulation for all three systems
# -----------------------------------------------------------------------
# initial condition
N = rep(N0,nr_spN)
P = rep(P0,nr_spP)
y0<-c(N=N,P=P)

res_corr<-lsoda(y=y0,times=times,func=Model_Prey_Pred,parms=para_corr)          # prey and predator with correlated niche axes
res_uncorr<-lsoda(y=y0,times=times,func=Model_Prey_Pred,parms=para_uncorr)      # prey and predator with uncorrelated niche axes
res_comp<-lsoda(y=y0[c(1:nr_spN)],times=times,func=Model_Comp,parms=para_corr)  # competition between species without predation
rm(para_corr,para_uncorr,K,sigma_spN,sigma_spP,c_,r,m_pred,e_pred,
   a_N,a_P_corr,a_P_uncorr,mu_N,mu_P,nu_N,N0,P0,N,P,y0)

# -----------------------------------------------------------------------
# Estimate variables - for one system at the time
# -----------------------------------------------------------------------
# cluster analysis
treshold <- 0.00001 # cluster quantification parameters

# format ode-results for further analysis
# NOTE: run the analysis for one of the systems at the time
(datForClust <- MakedatForClust(res_corr,muList,niche="corr",treshold=treshold))
(datForClust <- MakedatForClust(res_uncorr,muList,niche="uncorr",treshold=treshold))
(datForClust <- MakedatForClust(res_comp,muList,niche=NULL,treshold=treshold))
rm(res_corr,res_uncorr,res_comp)

datClust <- runKmeansGap(datForClust) # run cluster analysis for each community available in "datForClust" with pre-set function arguments
# Cluster algorithm originally from D'Andrea et al (2019) PLoS Comp. Biol. 15(1): e1006688. https://github.com/rafaeldandrea/Clustering-metric
# See end of script for running the KmeansGapMod() function directly
# NOTE: the clustering results can vary between runs because the function is based on an optimization algorithm and on generation of a random null model.
nicheplotClust(datClust=datClust,TL=2)
nicheplotClust(datClust=datClust,TL=1) # plot the clusters: TL = Number of Trophic Levels. (=2 for res_corr and res_uncorr, and =1 for res_1sp).

spp <- names(datClust) 
datRich <- sapply(spp, function(sp) {
  # calculate overall trait richness (OAtNr) and extinction (OAtExtProp)
  OAtNr <- sum(datForClust[[sp]][,"N"] > 0)
  nrStart <- get(paste0("nr_sp",if(sp == "M") "N" else sp))
  OAtExtProp <- (nrStart - OAtNr)/nrStart
  # calculate mean trait richness/cluster (BWtNr)
  IDuniq <- unique(datClust[[sp]]$clusterID[,"clustID"])
  BWtNr <- mean(sapply(IDuniq,function(id) {
    datID <- datClust[[sp]]$clusterID[which(datClust[[sp]]$clusterID[,"clustID"] == id),]
    traitNrClust <- sum(datID[,"abu"] > 0)
    return(traitNrClust)
  }))
  traitRich <- cbind(OAtNr=OAtNr,OAtExtProp=OAtExtProp,BWtNr=BWtNr)
  return(traitRich)
}, simplify = FALSE,USE.NAMES = TRUE)

datRes <- lapply(spp, function(sp) { 
  # gather results from datClust and datRich
  clust <- datClust[[sp]]$clustSum
  rich <- datRich[[sp]]
  dat <- cbind(clust,rich)
  colnames(dat) <- paste0(sp,".",colnames(dat))
  return(dat)
})


# Results: a vector of results for simulation i
###############################################
(SumRes <- c(if(length(spp)==1) para_sim_comp else para_sim,unlist(datRes)))
# The vector of results contain:
# 1) parameter value specific for the simulation
#     nSim      = ID of simulation, corresponding to the parameter values given in the object "parameters"
#     K         = Carrying capacity
#     c_        = Maximum attack rate
#     sigma_spN = Niche width of competing prey
#     sigma_spP = Niche width of predators
#     r         = Maximum per capita growth rate
#     m         = Mortality of predators 
#     e         = Conversion factor from one trophic level to the next
#     N0        = Abundance of prey at time step 0
#     P0        = Abundance of predators at time step 0
# 2) results. Depending on model ran above (corr, uncorr or comp), the vector will contain results for 
# N. = prey community with correlated niche axis, P. = predator community, and M. = prey community with correlated niche axis:
#   khat        = Number of clusters estimated for the observed community, from the KmeansGap algorithm
#   maxgap      = Gap statistic = maximum value of the gap index across all number of clusters tested in the observed community, from the KmeansGap algorithm
#   zscore      = Defined as (maxgap-mean(maxnullgap))/sd(maxnullgap), from the KmeansGap algorithm
#   pvalue      = Fraction of null communities with a higher gap statistic than the observed community, from the KmeansGap algorithm
#   OAtNr       = Overall number of species (richness)
#   OAtExtProp  = Overall proportion of extinction
#   BWtNr       = Mean number of species between clusters 




# -----------------------------------------------------------------------
# Results for all successful simulations and cluster identifications
# -----------------------------------------------------------------------
load("datCorr.R") # prey (N) and predators (P) with correlated niche axes for prey
load("datUncorr.R") # prey (N and M) and predators (P) with uncorrelated niche axes for prey
load("datComp.R") # competing community only (N)



# -----------------------------------------------------------------------
# Run the modified KmeansGapMod model directly
clust_Mod <- KmeansGapMod(dat=data,
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

# -----------------------------------------------------------------------