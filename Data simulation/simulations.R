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
    n.species=2, # The number of species
    n.sites=1000, # The number of sites
    n.visit=1 # The number of visits
    ),
fixed_effects = list(
  beta0 = c(1,0),#,0), # intercept for each species
  beta1=c(0,8),#,-2), #Covariate effect for each species,
  gamma0 = matrix(c(2, 1, 0,
                     1, 2, 0), nrow = 2,  
                   ncol = 3, byrow = TRUE),
  gamma1 = matrix(c(5, - 1,0 ,
                     -2, 4, 0), nrow = 2,  
                   ncol = 3, byrow = TRUE)
),
covariates = list(
  cov = runif(1000, -1,1), #NB: n.site = 400,
  cov_omega = runif(1000, -1,1)
)
)

#True values for the contrast of parameters for the coverage estimates
#in NIMBLE
true1 <- t(apply(input$fixed_effects$gamma0,1, 
      function(x){x[1:(length(x)-1)] - x[3]}))

true2 <- t(apply(input$fixed_effects$gamma1,1, 
        function(x){x[1:(length(x)-1)] - x[3]}))

## Generating data
genData <- function(input, type,seed=1){ # made all of the inputs null so can see everything feeds through
  #Retrieving Parameter values
  n.sites = input$constants$n.sites # number of sites
  n.visit=input$constants$n.visit #number of visits
  n.species=input$constants$n.species # Number of species
  beta0 <- input$fixed_effects$beta0 # intercept
  beta1 <- input$fixed_effects$beta1 # Covariate effect
  gamma0 <- input$fixed_effects$gamma0
  gamma1 <- input$fixed_effects$gamma1
  cov <- input$covariates$cov #covariate
  cov_omega <- input$covariates$cov_omega
  n.reported <- n.species +1
  n.validate <- floor(n.sites/3)
  #Check if parameter values are right
  if(length(n.sites) != 1) stop("Dimension of number of sites is wrong")
  if(length(n.visit) != 1) stop("Dimension of number of visits is wrong")
  if(length(n.species) != 1) stop("Dimension of number of species is wrong")
  if(length(beta0) != n.species) stop("Dimension of number of intercept is wrong")
  if(length(beta1) != n.species) stop("Dimension of number of covariate effects is wrong")
  if(length(cov) != n.sites) stop("Dimension of number of covariates is wrong")
  if(length(cov_omega) != n.sites) stop("Dimension of number of covariates is wrong")
  
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


#Verified species
library(boot)
library(extraDistr)
library(MASS)
for(site.tag in 1:n.sites){
  for(visit.tag in 1:n.visit){
    C[visit.tag, site.tag] <- rcat(1, prop[1:n.species,site.tag])
  }
}

alpha <- omega <- array(NA, dim = c(n.species, n.reported, n.sites))
if(type == "variable"){

for(site.tag in 1: n.sites){
  for(spe.tag in 1:n.species){
    for(report.tag in 1:(n.reported)){
      alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + gamma1[spe.tag, report.tag] * cov_omega[site.tag])
    }
  }
}
# Confusion Matrix for the Species
for(site.tag in 1: n.sites){
  for(spe.tag in 1:n.species){
    for(report.tag in 1:n.reported){
      omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
    }
  }
}
}


if(type == "constant"){
  for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      for(report.tag in 1:(n.reported)){
        alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag])
      }
    }
  }
  # Confusion Matrix for the Species
  for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      for(report.tag in 1:n.reported){
        omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
      }
    }
  }
}

for(visit.tag in 1:n.visit){
  for(site.tag in 1:n.sites){
    Y[visit.tag,site.tag] <- rcat(1,omega[C[visit.tag,site.tag],1:n.reported, site.tag])
  }
}


#Subsetting training and validation set
index_for_splitting <- sort(sample(1:n.sites, n.validate, replace = FALSE))
Validation_C <- C[, index_for_splitting]
Validation_Y <- Y[, index_for_splitting]
validation_indices_for_mismatch = index_for_splitting[which(Validation_C != Validation_Y)]
validation_indices_for_match = index_for_splitting[which(Validation_C == Validation_Y)]
validation_mismatchC = C[, validation_indices_for_mismatch]
validation_mismatchY = Y[, validation_indices_for_mismatch]
validation_matchC = C[, validation_indices_for_match]
validation_matchY = Y[, validation_indices_for_match]
C[1:n.visit, index_for_splitting] <- NA



  dataout <- list(C=C, #Verified data
                  Y=Y, # Reported data
                  Lam=lambda, #Mean intensity
                  omega=omega, # Confusion matrix
                  cov=cov, # covariate forecological process
                  cov_omega = cov_omega, #covariate for the classification process
                  prop=prop, #proportion of intensity,
                  true_values = c(true1, true2), #estimates of contrasts of classification process pars
                  Validation_C = Validation_C, #verified data for validation
                  Validation_Y = Validation_Y, #reported data for validation
                  validation_indices_for_mismatch = validation_indices_for_mismatch, # IDs for the data with mismatch 
                  index_for_splitting = index_for_splitting, #IDs for the data for validation
                  validation_mismatchC = validation_mismatchC, #verified data with mismatched classifications
                  validation_mismatchY = validation_mismatchY, #reported data with mismatched classifications
                  validation_indices_for_match = validation_indices_for_match, #IDs for data for correctly matched
                  validation_matchC = validation_matchC, #verified data of correctly matched observations
                  validation_matchY = validation_matchY #reported data of correctly matched observations
                  )
  return(dataout)
  
}

#Simulate 100 replicates by changing the seed value for each replicate
nreplicates <- 100

#Heterogeneous classification process
sim1 <- pblapply(1:nreplicates, function(x){
  data <- genData(input, seed = x, type = "variable")
}, cl=1)

#Homogeneous classification process
sim2 <- pblapply(1:nreplicates, function(x){
  data <- genData(input, seed = x, type = "constant")
}, cl=1)

sim <- c(sim1, sim2)

save(sim, file="simulated data1.RData")

