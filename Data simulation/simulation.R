#Load packages needed for this package to run
library(spatstat)
library(pbapply)
library(boot)
library(extraDistr)
library(MASS)

#Set seed for the simulation
set.seed(1994)

# Setting parameter values needed for the model
input <- list(
  constants=list(
    n.species=3, # The number of species
    n.sites=400, # The number of sites
    n.visit=1, # The number of visits
    p.tag= c(0.7,0.8,0.90) # The probability of correct identification for the confusion matrix
),
fixed_effects = list(
  beta0 = c(-3,1,2), # intercept for each species
  beta1=c(1,0,-1) #Covariate effect for each species
),
covariates = list(
  cov = runif(400, -2,2) #NB: n.site = 400
)
)

## Generating data
genData <- function(input, seed=1){ # made all of the inputs null so can see everything feeds through
  #Retrieving Parameter values
  n.sites = input$constants$n.sites # number of sites
  n.visit=input$constants$n.visit #number of visits
  p.tag = input$constants$p.tag # probability of correct classification
  n.species=input$constants$n.species # Number of species
  beta0 <- input$fixed_effects$beta0 # intercept
  beta1 <- input$fixed_effects$beta1 # Covariate effect
  cov <- input$covariates$cov #covariate
  
  #Check if parameter values are right
  if(length(n.sites) != 1) stop("Dimension of number of sites is wrong")
  if(length(n.visit) != 1) stop("Dimension of number of visits is wrong")
  if(length(n.species) != 1) stop("Dimension of number of species is wrong")
  if(length(beta0) != n.species) stop("Dimension of number of intercept is wrong")
  if(length(beta1) != n.species) stop("Dimension of number of covariate effects is wrong")
  if(length(cov) != n.sites) stop("Dimension of number of covariates is wrong")
  
  #Estimating mean intensity
  lambda <- matrix(NA, nrow=n.species, ncol=n.sites)
  for(spe.tag in 1:n.species){
    for(site.tag in 1:n.sites){
      lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag]*cov[site.tag])
    }
  }

#Proportion of each species as defined by Equation (2) in the paper
 prop <-  matrix(NA, ncol=n.sites ,nrow=n.species) 
for(site.tag in 1:n.sites){
  prop[,site.tag]<- proportions(lambda[1:n.species, site.tag])
}

#Check if the proportions are correct and sum up to 1
for(site.tag in n.sites){
  for(spe.tag in 1:n.species){
    if(any(prop[spe.tag, site.tag] <0) | any(prop[spe.tag, site.tag] >1) ) stop("Entry in Probability wrong")
    if(round(sum(prop[1:n.species,site.tag ]),0) != 1) stop("Sum of probabilities is not equal to 1")
  }
}

z <- array(NA, dim=c(n.species, n.visit, n.sites)) #Detect/Non-Detection
C <- array(NA, dim=c(n.visit,n.sites)) # Obs detect/non-detect
Y <- numb <- numb_with_NA <- array(NA, dim=c(n.visit, n.sites)) # Verified species with Machine Learning
omega <- array(NA, dim= c(n.species, n.species+1)) # Classification Confusion matrix

for(site.tag in 1:n.sites){
  #   for(visit.tag in 1:n.visit){
  omega[1,2] <- 0.05
  omega[1,3] <- 0.13
  omega[2,1] <- 0.10
  omega[2,3] <- 0.05
  omega[3,1] <- 0
  omega[3,2] <- 0
  for(spe.tag in 1:n.species){
    omega[spe.tag, spe.tag] <- p.tag[spe.tag]
    omega[spe.tag, n.species+1] <- 1 - sum(omega[spe.tag, 1:n.species])  
    #}
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

#Verified species
library(boot)
library(extraDistr)
library(MASS)
for(site.tag in 1:n.sites){
  for(visit.tag in 1:n.visit){
    C[visit.tag, site.tag] <- rcat(1, prop[1:n.species,site.tag])
  }
}

# Simulating the reported data
for(visit.tag in 1:n.visit){
  for(site.tag in 1:n.sites){
    Y[visit.tag,site.tag] <- rcat(1, prob = omega[C[visit.tag,site.tag], 1:(n.species+1)])
  }
}

  dataout <- list(C=C, #Verified data
                  Y=Y, # Reported data
                  Lam=lambda, #Mean intensity
                  omega=omega, # Confusion matrix
                  cov=cov, # covariate
                  prop=prop #proportion of intensity
                  )
  return(dataout)
  
}

#Simulate 1000 replicates by changing the seed value for each replicate
nreplicates <- 1000
sim <- pblapply(1:nreplicates, function(x){
  data <- genData(input, seed = x)
}, cl=1)


save(sim, file="simulated data1.RData")
