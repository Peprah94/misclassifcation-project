#This code is used to estimate the parameters for the model from the simulated data

#Load the simulated data from the "simulation.R" code
load("simulated data1.RData")

# Load packages needed to run the code
library(pbapply)
library(nimble)
library(MCMCglmm)
library(coda)
library(mcmcplots)
library(MCMCpack)
library(boot)
library(extraDistr)
library(MASS)
library(spdep)
library(parallel)
require(ggmcmc)

# Function that estimates the parameters
#Needed to run the code in parallel
run_simulations <- function(simulations_all){
  start_time <- Sys.time()
  
  # packages copied within the code for the parallelisation
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(mcmcplots)
  library(MCMCpack)
  library(boot)
  library(extraDistr)
  library(MASS)
  library(spdep)
  require(ggmcmc)

code <- nimbleCode({
  # Priors for beta
  for(spe.tag in 1:n.species){
 beta0[spe.tag] ~ dnorm(0,sd=1)
    beta1[spe.tag] ~ dnorm(0, sd=1)
  }

  # Linear Predictor
  for(spe.tag in 1:n.species){
    for(site.tag in 1:n.sites){
      linpred[spe.tag,site.tag] <-  beta0[spe.tag] + beta1[spe.tag]*cov[site.tag] 
     lambda[spe.tag,site.tag] <- exp(linpred[spe.tag,site.tag])
    }
  }
  
  # Proportion for the multinomial distribution
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
    prop[spe.tag,site.tag] <- (lambda[spe.tag, site.tag])/sum(lambda[1:n.species, site.tag])
    }
  }
  
  
  #Prior for the Confusion Matrix
  for(spe.tag in 1:n.species){
    for(i in 1:(n.species+1)){
      alpha[spe.tag, i] ~ dexp(1)
    }
  }
  # Confusion Matrix for the Species
  for(spe.tag in 1:n.species){
    omega[spe.tag, 1:(n.species+1)] ~ ddirch(alpha[spe.tag,1:(n.species+1)])
  }
  
  # Verified data
  for(site.tag in 1:n.sites){
    for(visit.tag in 1:n.visit){
      C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
    }
  }

  # Reported observations
  for(visit.tag in 1:n.visit){
    for(site.tag in 1:n.sites){
      Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:(n.species+1)])
    }
  }
  #estimating contrasts
      z021 <- beta0[2]- beta0[1]
    z031 <- beta0[3] - beta0[1]
    z121 <- beta1[2]- beta1[1]
    z131 <- beta1[3] - beta1[1]
})

# Retrieving data from the simulated data
C <- simulations_all$C # Unverified citizen science observation
Y <- simulations_all$Y #Verified Citizen science observation
cov <- simulations_all$cov #Covariates
data <- list(C,Y, cov)

#Constants for the model
dim_data <- dim(data[[1]])
n.cov=ncol(cov)
const <- list(n.sites=dim_data[2],
              n.species=max(c(data[[1]])),
              n.visit=(dim_data[1]), 
              n.cov=n.cov
              )

# Data for the Model
idm_data <- list(C = data[[1]], 
                 Y = data[[2]], 
                 cov=data[[3]])

# Initials for the model
alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.species+1))
for(spe.tag in 1:const$n.species){
  for(k in 1:(const$n.species+1)){
    alpha[spe.tag,k] <- 1
  }
}
omega <- matrix(NA, nrow=const$n.species,
                ncol = (const$n.species)+1) 
for(spe.tag in 1:(const$n.species)){
  omega[spe.tag,]<- rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.species+1)])
}

idm_inits <- function(){list(beta1 =rep(1, const$n.species),
                             omega=omega,
                             beta0 = rep(1, const$n.species),
                             alpha = alpha)
}
initsList <- idm_inits()

#Putting all together for the creating compilation
modelInfo <- list(
  code = code,
  constants = const,
  data = idm_data,
  inits = initsList
)

#Create the model in nimble
mwtc <- nimbleModel(code,data = idm_data, 
                    constants = const, 
                    inits = initsList,
                    dimensions =list(lambda = c(const$n.species, const$n.sites), 
                                     calculate = FALSE))

# Create the model in C
Cmwtc <- compileNimble(mwtc, 
                       showCompilerOutput = FALSE)

  mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0","beta1","omega", "z021", "z031", "z121", "z131"))

Rmcmc <- buildMCMC(mcmcconf, 
                   enableWAIC =FALSE)

# Compile 
cmcmc <- compileNimble(Rmcmc, 
                       project = Cmwtc,
                       resetFunctions = TRUE)

# Run the MCMC
mcmc.out <- runMCMC(cmcmc, 
                    niter = 100000,
                    nchains = 3,
                    nburnin = 50000,
                    inits = initsList,
                    thin =10, 
                    setSeed = TRUE, 
                    samples=TRUE, 
                    samplesAsCodaMCMC = TRUE, 
                    summary = TRUE, 
                    WAIC = FALSE)

# Return the summary of MCMC results
output <- mcmc.out$summary$all.chains
end.time <- Sys.time()
time_taken <- as.numeric(round(end.time-start_time,2))

# Preparing the results for the extraction of Rhat
mcmclist <- ggs(mcmc.out$samples)

# Subset the omega and beta parameters
aa <- mcmclist%>%
  filter(grepl('omega|beta',Parameter))

subset_parameters<- unique(aa$Parameter)

# Rhat values
subset_Rhat <- aa%>%
  ggs_Rhat()
Rhat_data <- subset_Rhat$data[,5]
all_rhat <- all(Rhat_data < 1.05) #Rhat is less than 1.05
N_over_rhat <-length(which(Rhat_data > 1.05))/length(subset_parameters) #Rhats over 1.04

  #estimate coverage
coverage <- function(output){
    ifelse(output[4] <= output[6] & output[5] >= output[6],1,0)
  }
  true_values <- c(0.7, 0.1, 0, 0.05,0.8, 0, 0.13, 0.05, 0.9, 0.12,0.05,0.1,
                   2,1,-2,-4)
  df <- cbind(output[7:22,], true_values)
  coveragevalues = apply(df,1, coverage)
  returnedcoverages <- c(coveragevalues,N_over_rhat=N_over_rhat, all_rhat = all_rhat)
  
  #return results
  beta0 <- output[(1:const$n.species),2]
  beta1 <- output[((const$n.species+1):(const$n.species*2)),2]
  #lambda<- matrix(output[((const$n.species*2)+1) : (((const$n.species*2))+ (const$n.species*(const$n.sites))),2], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  omega <- matrix(output[((const$n.species*2)+ 1): ((((const$n.species*2)))+(const$n.species*(const$n.species+1))),2], nrow=const$n.species, ncol=(const$n.species+1), byrow=FALSE)
  contrasts <- tail(output[,2],4)
  return(list(omega=omega, beta0=beta0, beta1 = beta1,contrasts=contrasts, returnedcovarages=returnedcoverages))
}

#Number of clusters
cl <- makeCluster(5)
setDefaultCluster(cl)
rep_estimates <- pblapply(sim, run_simulations, cl=cl)
stopCluster(cl)

#Save results
save(rep_estimates, file="estimateddata_new.RData")

