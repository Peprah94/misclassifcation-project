# Set-up for PACKAGES
if (!require(nimble)) install.packages('nimble')
if (!require(MCMCglmm)) install.packages('MCMCglmm')
if (!require(coda)) install.packages('coda')
if (!require(mcmcplots)) install.packages('mcmcplots')
if (!require(MCMCpack)) install.packages('MCMCpack')
if (!require(boot)) install.packages('boot')
if (!require(extraDistr)) install.packages('extraDistr')


#source("functions.R")
library(pbapply)

#Function that simulates data
simrep <- function(input,seeds){
  n.species=input$n.species #Number of species
  n.sites=input$n.sites  #Number of sites
  n.visit = input$n.visit #Number of visits
  n.cov=input$n.cov #Number of covariates
  p.tag = input$p.tag  #Detection probability
  beta=input$beta  # Coefficients of the model
  #
  set.seed(seeds)
  
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(mcmcplots)
  library(MCMCpack)
  library(boot)
  library(extraDistr)
  library(MASS)
  library(spdep)

  
  source("maternfnx.R")
  #Dimensions of data
  cov <- array(NA, dim=c(n.sites, n.cov)) #Covariates
  lambda <- prop <- array(NA, dim=c(n.species,n.sites)) #Occupancy prob
  z <- array(NA, dim=c(n.species, n.visit, n.sites)) #Detect/Non-Detection
  C <- array(NA, dim=c(n.species, n.visit,n.sites)) # Obs detect/non-detect
  Y <- numb <- array(NA, dim=c(n.visit, n.sites)) # Verified species with Machine Learning
  omega <- array(NA, dim= c(n.species, n.species+1)) # Classification Confusion matrix
  
  #check dimensions of beta
  if(dim(beta)[1]!= n.species | dim(beta)[2]!= (n.cov+1)) stop("Dimensions of beta are wrong")
  
  # Random effects for sites
  #sites <- seq(1,n.sites,1)
 # alpha_sites <- rnorm(n.sites, mean = 0, sd= 1)
  x0 = runif(n.sites, 0, 3)
  y0 = runif(n.sites, 0, 3)
  locs <- cbind(x0,y0)
  V=cov.matern(as.matrix(dist(locs)),nu=1/2, alpha=7*sqrt(3))
alpha_sites=mvrnorm(mu=rep(0, n.sites), Sigma=V)

  #Simulation of covariates
cov <- matrix(runif(n.sites*n.cov, -1,1), nrow = n.sites, ncol = n.cov)
  
  # Linear Predictor the P
  for(spe.tag in 1:n.species){
    for(site.tag in 1:n.sites){
      #lambda[spe.tag,site.tag] <- exp(beta[spe.tag,1] + alpha_sites[site.tag])
  lambda[spe.tag,site.tag] <- exp(beta[spe.tag,1] + beta[spe.tag,2]*cov[site.tag,1]+ alpha_sites[site.tag])
  }
  }

  # Detection/ Non-detection
  for(site.tag in 1:n.sites){
      prop[ ,site.tag]<- proportions(lambda[1:n.species, site.tag])
  }
  
#Check if the probabilities are correct and sum up to 1
    for(site.tag in n.sites){
      for(spe.tag in 1:n.species){
        if(any(prop[spe.tag, site.tag] <0) | any(prop[spe.tag, site.tag] >1) ) stop("Entry in Probability wrong")
        if(round(sum(prop[1:n.species,site.tag ]),0) != 1) stop("Sum of probabilities is not equal to 1")
      }
  }
  
  # Confusion Matrix
  for(site.tag in 1:n.sites){
    for(visit.tag in 1:n.visit){
      omega[1,2] <- 0.05
      omega[1,3] <- 0.13
      omega[2,1] <- 0.15
      omega[2,3] <- 0.1
      omega[3,1] <- 0.09
      omega[3,2] <- 0.18
      for(spe.tag in 1:n.species){
        omega[spe.tag, spe.tag] <- p.tag 
        omega[spe.tag, n.species+1] <- 1 - sum(omega[spe.tag, 1:n.species])  
      }
    }
  }
  
  # Check if the rowSums are 1
  for(spe.tag in 1:n.species){
    for(visit.tag in 1:n.visit){
      for(site.tag in n.sites){
        if(sum(omega[spe.tag,]) != 1) stop("Sum of confusion matrix is not equal to 1")
        if(any(omega[spe.tag,] <0) | any(omega[spe.tag,] >1) ) stop("Entry with Confusion matrix")
      }
    }
  }
  
  # Observed detections/non-detections
  for(site.tag in 1:n.sites){
    for(visit.tag in 1:n.visit){
      C[,visit.tag, site.tag] <- rmultinom(1,1, prop[1:n.species,site.tag])
    }
  }
  

# Which species were present
  for(visit.tag in 1:n.visit){
    for(site.tag in 1:n.sites){
      numb[visit.tag,site.tag] <-which(C[1:n.species,visit.tag,site.tag]==1)
    }
  }
  
  # Simulating the verified data
  for(visit.tag in 1:n.visit){
    for(site.tag in 1:n.sites){
      Y[visit.tag,site.tag] <- rcat(1, prob = omega[numb[visit.tag,site.tag],1:(n.species+1)])
    }
  }
  
neigh <- dnearneigh(locs, d1=0, d2=sqrt(1)*1000+1)
winnb <- nb2WB(neigh)

  # Returning the needed information  
  return(list(C=C,
              Y=Y, 
              omega=omega, 
              cov=cov, 
              prop=prop, 
              locs=locs, 
              winnb=winnb))
}


#Input
input <- list(
  n.species=3,
  n.sites=200,
  n.visit=1, 
  n.cov=1,
  beta = matrix(c(2,-1,
                  -1,3,
                  1,2),
                nrow = 3,
                ncol = 2),
  p.tag= 0.7
)


#Simulation of the data
#sim <- simrep(input, 1)
nreplicates <- 30
sim <- pblapply(1:nreplicates, function(x){
  data <- simrep(input, seeds = x)
}, cl=3)


save(sim, file="simulated data.RData")


