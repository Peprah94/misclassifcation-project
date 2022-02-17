#Load the data from the data manipulation
load("data_new.RData")

# Load packages
library(parallel)
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
library(ggmcmc)

# Function for running the data analysis
run_simulations <- function(simulations_all){
  #Load packages
  library(parallel)
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
  library(ggmcmc)
  
  code <- nimbleCode({
    # Priors for beta
    for(spe.tag in 1:n.species){
      beta0[spe.tag]~dnorm(0, sd=1)
      beta1[spe.tag] ~ dnorm(0, sd=1)
      beta2[spe.tag] ~ dnorm(0, sd=1)
      beta3[spe.tag] ~ dnorm(0, sd=1)
    }
    
    # Linear Predictor
    for(spe.tag in 1:n.species){
      for(site.tag in 1:n.sites){
        linpred[spe.tag,site.tag] <-  beta0[spe.tag] + beta1[spe.tag]*temp[site.tag]+ beta2[spe.tag]*altitudes[site.tag]+ beta3[spe.tag]*agree[site.tag]
        log(lambda[spe.tag,site.tag]) <- linpred[spe.tag,site.tag]
      }
    }
    
    # Proportion for the multinomial distribution
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
      }
    }
    
    
    #Prior for the Confusion Matrix
    for(spe.tag in 1:n.species){
      for(i in 1:n.reported){
        alpha[spe.tag, i] ~ dexp(1)
      }
    }
    # Confusion Matrix for the Species
    for(spe.tag in 1:n.species){
      omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
    }
    
    # Verified observations 
    for(site.tag in 1:n.sites){
      C[site.tag] ~ dcat(prop[1:n.species,site.tag])
    }
    
    # Verified observations from ML classifications
    for(site.tag in 1:n.sites){
      Y[site.tag] ~ dcat(omega[C[site.tag],1:n.reported])
    }
  })
  
  # Retrieving data from the simulations
  C <- (simulations_all$C) # Unverified citizen science observation
  Y <- (simulations_all$Y) #Verified Citizen science observation
  agree <- (simulations_all$agree)
  disagree <- simulations_all$disagree#Covariates
  temp <- simulations_all$temp #temperature
  altitudes <- simulations_all$altitude
  
 
  #Setting parameters
  data <- list(C,Y, temp, agree, disagree, altitudes)
  n.species <- length(unique(c(data[[1]])))
  n.reported <- length(unique(c(data[[2]])))
  n.sites <- length(data[[1]])

  
  # Constants for the MCMC model
  const <- list(n.sites=n.sites,
                n.species=n.species,
                n.reported=n.reported)
  
  alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.reported))
  for(spe.tag in 1:const$n.species){
    for(k in 1:(const$n.reported)){
      alpha[spe.tag,k] <- 1
    }
  }

  
  # Data for the Model
  idm_data <- list(C = data[[1]], 
                   Y = data[[2]], 
                   temp=data[[3]],
                   agree=data[[4]],
                   disagree=data[[5]],
                   altitudes = c((data[[6]])))
  
  # Initials for the model
  omega <- matrix(NA, nrow=const$n.species,
                  ncol = (const$n.reported)) 
  for(spe.tag in 1:(const$n.species)){
    omega[spe.tag,]<- rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.reported)])
  }
  
  idm_inits <- function(){list(beta1 =rep(1, const$n.species),
                               beta2 =rep(1, const$n.species),
                               beta3 =rep(1, const$n.species),
                               beta0 =rep(1, const$n.species),
                               omega=omega,
                               mean.lam = rep(1,const$n.species),
                               tau= 2,
                               alpha = alpha,
                               alpha_sites = rep(0, const$n.sites)
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
  
  mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0","beta1","beta2", "beta3","tau","omega", "prop","lambda"))
  Rmcmc <- buildMCMC(mcmcconf, 
                     enableWAIC =FALSE)
  
  # Compile 
  cmcmc <- compileNimble(Rmcmc, 
                         project = Cmwtc,
                         resetFunctions = TRUE)
  
  # Run the MCMC
  mcmc.out <- runMCMC(cmcmc, 
                      niter = 300000,
                      nchains = 3,
                      nburnin = 200000,
                      inits = initsList,
                      thin =50, 
                      setSeed = TRUE, 
                      samples=TRUE, 
                      samplesAsCodaMCMC = TRUE, 
                      summary = TRUE, 
                      WAIC = FALSE)
  
  mcmclist <- ggs(mcmc.out$samples)
  ggs_density(mcmclist, family = "omega")
  ggs_density(mcmclist, family = "beta")
  ggs_traceplot(mcmclist, family = "beta")
  ggs_traceplot(mcmclist, family = "omega")
  subset_parameters <- c("beta0[1]", "beta1[1]", "beta2[1]", "beta3[1]",
                         "beta0[2]", "beta1[2]", "beta2[2]", "beta3[2]",
                         "beta0[3]", "beta1[3]", "beta2[3]", "beta3[3]")
aa <-mcmclist%>%
    filter(Parameter %in% subset_parameters)%>%
  ggs_Rhat()
 Rhat_data <- aa$data[,5]
all_rhat <- all(Rhat_data < 1.05) #Rhat is less than 1.05
N_over_rhat <-length(which(Rhat_data > 1.05))/length(subset_parameters) #Rhats over 1.05
  
  output <- mcmc.out$summary$all.chains
  beta0 <- output[(1:const$n.species),1]
  beta1 <- output[((const$n.species+1):(const$n.species*2)),1]
  beta2 <- output[(((const$n.species*2)+1):(const$n.species*3)),1]
  beta3 <- output[(((const$n.species*3)+1):(const$n.species*4)),1]
  lambda<- matrix(output[((const$n.species*4)+1) : (((const$n.species*4))+ (const$n.species*(const$n.sites))),1], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  omega <- matrix(output[((((const$n.species*4))+ (const$n.species*(const$n.sites)))+1): ((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))),1], nrow=const$n.species, ncol=(const$n.reported), byrow=FALSE)
  prop <-matrix(output[((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))+1):((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))+(const$n.species*const$n.sites)),1],nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  tau <- tail(output[,1])[6]
  return(list(omega=omega, beta0=beta0, beta1 = beta1,beta2=beta2,beta3=beta3, prop=prop, tau=tau,lambda=lambda, mcmc.out=mcmc.out,all_rhat = all_rhat,N_over_rhat = N_over_rhat))
}

cl <- makeCluster(6)
setDefaultCluster(cl)
rep_estimates <- pblapply(data_for_nimble,run_simulations, cl=cl)
stopCluster(cl)
save(rep_estimates, file="data_results.RData")

