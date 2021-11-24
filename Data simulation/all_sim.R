load("simulated data1.RData")
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


#source("maternfnx.R")

#index <- function(y){
 # x <- nimC(y)
#  return(which(x==1))
#}
#Rindex <- nimbleRcall(function(y=double(1)){}, Rfun="index", returnType = double(0))

# Nimble function for estimating the proportions
#nim_proportion <- function(y){
  #x <- nimC(y)
  #z <- (x)
  #ret <- proportions(z)
 # return(ret)
#}
#nimbleProportion <- nimbleRcall(function(y=double(1)){}, Rfun="nim_proportion", returnType = double(1))

run_simulations <- function(simulations_all){
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(mcmcplots)
  library(MCMCpack)
  library(boot)
  library(extraDistr)
  library(MASS)
  library(spdep)


  #source("maternfnx.R")
  
  #index <- function(y){
   # x <- nimC(y)
   # return(which(x==1))
 # }
 # Rindex <- nimbleRcall(function(y=double(1)){}, Rfun="index", returnType = double(0))
  
  # Nimble function for estimating the proportions
 # nim_proportion <- function(y){
  #  x <- nimC(y)
   #z <- (x)
    #ret <- proportions(z)
    #return(ret)
  #}
  #nimbleProportion <- nimbleRcall(function(y=double(1)){}, Rfun="nim_proportion", returnType = double(1))

#The code for the MCMC
code <- nimbleCode({
  
  
  # Priors for beta
  for(spe.tag in 1:n.species){
    #for(cov.tag in 1:(n.cov+1)){
    #beta[spe.tag,cov.tag] ~ dnorm(0, sd=100)
    #beta[spe.tag,1] ~ dnorm(0, sd=100)
    #mean.lam[spe.tag]~dunif(0,20)
    beta0[spe.tag] ~ dnorm(0,sd=1)
   # beta0[spe.tag]~dnorm(0, sd=10)
    beta1[spe.tag] ~ dnorm(0, sd=1)
    #beta2[spe.tag] ~ dnorm(0, sd=1)
    #}
  }
  
  
  #Prior for spatial field
  #tau ~dgamma(1,1) #Chi-square
  #v.eta ~ T(dnorm(0, 0.1), 0, ) # I(0,)
  #tau <- 1/v.eta
  
  
  # Random Effects for the model
  # sig_sites ~ dgamma(1,1)
  # mu_sites ~ dnorm(0, 1)
  #for(site.tag in 1:n.sites){
  #alpha_sites[site.tag] ~ dnorm(mu_sites, sig_sites)
  #alpha_sites[1:n.sites] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.sites], tau, zero_mean = 1)  
  #}
  
  # Linear Predictor
  for(spe.tag in 1:n.species){
    for(site.tag in 1:n.sites){
      #log(lambda[spe.tag,site.tag]) <- beta[spe.tag,1] +  alpha_sites[site.tag]
      linpred[spe.tag,site.tag] <-  beta0[spe.tag] + beta1[spe.tag]*cov[1,site.tag] #+ beta2[spe.tag]*cov[2,site.tag]
      #lim.linpred[spe.tag,site.tag] <- min(1000, max(-1000, linpred[spe.tag,site.tag])) #stabilize the log
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
  
  # Citizen Science observations
  for(site.tag in 1:n.sites){
    for(visit.tag in 1:n.visit){
      C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
    }
  }
  
  # Subsetting index of which species are present
#for(visit.tag in 1:n.visit){
 # for(site.tag in 1:n.sites){
 #   numb[visit.tag,site.tag] <- C[visit.tag, site.tag]
 # }
 #}
  
  # Verified observations from ML classifications
  for(visit.tag in 1:n.visit){
    for(site.tag in 1:n.sites){
      Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:(n.species+1)])
    }
  }
  
})

# Retrieving data from the simulations
C <- simulations_all$C # Unverified citizen science observation
Y <- simulations_all$Y #Verified Citizen science observation
cov <- simulations_all$cov #Covariates

#Preparing the spatial component of the data
#locs <- cbind(simulations_all$locs$Var1, simulations_all$locs$Var2)
#neigh <- dnearneigh(locs, d1=0, d2=sqrt(0.5)+1)
#winnb <- nb2WB(neigh)
#winnb <- simulations_all$winnb
#L <- length(winnb$adj)
#str(winnb)
#Setting parameters
data <- list(C,Y, cov)
dim_data <- dim(data[[1]])
n.cov=ncol(cov)

n.species <- max(c(C))
#alpha <- rep(1, (n.species+1))

# Constants for the MCMC model
const <- list(n.sites=dim_data[2],
              n.species=max(c(data[[1]])),
              #n.sites = dim_data[3], 
              #n.species = (dim_data[1]), 
              n.visit=(dim_data[1]), 
              # n.visit=1, #use only for data
              n.cov=n.cov
              #L=L, 
              #adj=winnb$adj, 
              #num=winnb$num, 
              #weights=winnb$weights
              )
alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.species+1))
for(spe.tag in 1:const$n.species){
  for(k in 1:(const$n.species+1)){
    alpha[spe.tag,k] <- 1
 }
}
#alpha <- rep(1, (const$n.species+1))

# Data for the Model
idm_data <- list(C = data[[1]], 
                 Y = data[[2]], 
                 cov=data[[3]])

# Initials for the model
omega <- matrix(NA, nrow=const$n.species,
                ncol = (const$n.species)+1) 
for(spe.tag in 1:(const$n.species)){
  omega[spe.tag,]<- rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.species+1)])
}

idm_inits <- function(){list(beta1 =rep(1, const$n.species),
                             #beta2 = rep(1,const$n.species),
                             omega=omega,
                             beta0 = rep(1, const$n.species),
                             #tau= 2,
                             alpha = alpha
                             #alpha_sites = rep(0, const$n.sites )
)
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


mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0","beta1","omega", "prop","lambda"))

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

output <- mcmc.out$summary$all.chains
beta0 <- output[(1:const$n.species),2]
beta1 <- output[((const$n.species+1):(const$n.species*2)),2]
lambda<- matrix(output[((const$n.species*2)+1) : (((const$n.species*2))+ (const$n.species*(const$n.sites))),2], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
omega <- matrix(output[((((const$n.species*2))+ (const$n.species*(const$n.sites)))+1): ((((const$n.species*2))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.species+1))),2], nrow=const$n.species, ncol=(const$n.species+1), byrow=FALSE)
#prop <- matrix(output[(((((const$n.species*2))+ (const$n.species*(const$n.sites)))+(const$n.species*(n.species+1)))) :((((const$n.species*2))+ (const$n.species*(const$n.sites)))+(const$n.species*(n.species+1)))+ (const$n.species*(const$n.sites)),1], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
prop <-matrix(tail(output[,2],(const$n.species*(const$n.sites))),nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
return(list(omega=omega, beta0=beta0, beta1 = beta1, prop=prop, lambda=lambda))
}


cl <- makeCluster(5)
setDefaultCluster(cl)
rep_estimates <- pblapply(sim, run_simulations, cl=cl)
#rep_estimates <- pblapply(sim[[1]], function(x){
#  est <- run_simulations(x)
#}, cl=1)
#rep_estimates <- run_simulations(sim[[1]])

save(rep_estimates, file="estimateddata_new.RData")

