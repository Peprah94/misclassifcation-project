load("data_new.RData")
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

#handlers(global = TRUE)


#source("maternfnx.R")

index <- function(y){
  x <- nimC(y)
  return(which(x==1))
}
Rindex <- nimbleRcall(function(y=double(1)){}, Rfun="index", returnType = double(0))

#Nimble function for estimating the proportions
nim_proportion <- function(y){
  x <- nimC(y)
  z <- (x)
  ret <- proportions(z)
  return(ret)
}

nimbleProportion <- nimbleRcall(function(y=double(1)){}, Rfun="nim_proportion", returnType = double(1))

run_simulations <- function(simulations_all){
  
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
  
  #source("maternfnx.R")
  
  # index <- function(y){
  #  x <- nimC(y)
  #  return(which(x==1))
  # }
  # Rindex <- nimbleRcall(function(y=double(1)){}, Rfun="index", returnType = double(0))
  
  # Nimble function for estimating the proportions
  # nim_proportion <- function(y){
  #x <- nimC(y)
  #  z <- (x)
  # ret <- proportions(z)
  # return(ret)
  # }
  #nimbleProportion <- nimbleRcall(function(y=double(1)){}, Rfun="nim_proportion", returnType = double(1))
  
  #The code for the MCMC
  code <- nimbleCode({
    
    
    # Priors for beta
    for(spe.tag in 1:n.species){
      #for(cov.tag in 1:(n.cov+1)){
      #beta[spe.tag,cov.tag] ~ dnorm(0, sd=100)
      #beta[spe.tag,1] ~ dnorm(0, sd=100)
      #mean.lam[spe.tag]~dunif(0,20)
      #beta0[spe.tag] <- log(mean.lam[spe.tag]) 
      beta0[spe.tag]~dnorm(0, sd=1)
      beta1[spe.tag] ~ dnorm(0, sd=1)
      beta2[spe.tag] ~ dnorm(0, sd=1)
      beta3[spe.tag] ~ dnorm(0, sd=1)
      #beta4[spe.tag] ~ dnorm(0, sd=1)
      #}
    }
    
    
    #Prior for spatial field
    tau ~dgamma(1,1) #Chi-square
    
    
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
       #linpred[spe.tag,site.tag] <-  beta0[spe.tag] + beta1[spe.tag]*temp[site.tag]+ beta2[spe.tag]*agree[site.tag]+beta3[spe.tag]*disagree[site.tag]+ beta4[spe.tag]*altitudes[site.tag]#+alpha_sites[site.tag]
        #linpred[spe.tag,site.tag] <-  beta0[spe.tag] + beta1[spe.tag]*temp[site.tag]+ beta2[spe.tag]*altitudes[site.tag]#lim.linpred[spe.tag,site.tag] <- min(1000, max(-1000, linpred[spe.tag,site.tag]))
        linpred[spe.tag,site.tag] <-  beta0[spe.tag] + beta1[spe.tag]*temp[site.tag]+ beta2[spe.tag]*altitudes[site.tag]+ beta3[spe.tag]*agree[site.tag]#lim.linpred[spe.tag,site.tag] <- min(1000, max(-1000, linpred[spe.tag,site.tag])) #stabilize the log
        log(lambda[spe.tag,site.tag]) <- linpred[spe.tag,site.tag]
      }
    }
    
    # Proportion for the multinomial distribution
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        #prop[1:n.species,site.tag] <- nimbleProportion(lambda[1:n.species, site.tag]) 
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
    
    # Citizen Science observations
    for(site.tag in 1:n.sites){
      #for(visit.tag in 1:n.visit){
      C[site.tag] ~ dcat(prop[1:n.species,site.tag])
      #}
    }
    
    # Subsetting index of which species are present
    #for(visit.tag in 1:n.visit){
    #  for(site.tag in 1:n.sites){
    #    numb[visit.tag,site.tag] <- Rindex(C[1:n.species,visit.tag, site.tag])
    #  }
    #}
    
    # Verified observations from ML classifications
    #for(visit.tag in 1:n.visit){
    for(site.tag in 1:n.sites){
      Y[site.tag] ~ dcat(omega[C[site.tag],1:n.reported])
    }
    # }
    
  })
  
  # Retrieving data from the simulations
  C <- (simulations_all$C) # Unverified citizen science observation
  Y <- (simulations_all$Y) #Verified Citizen science observation
  #cov <- c(t(simulations_all$cov)) 
  agree <- (simulations_all$agree)
  disagree <- simulations_all$disagree#Covariates
  temp <- simulations_all$temp
  altitudes <- simulations_all$altitude
  
  #Preparing the spatial component of the data
  #locs <- cbind(simulations_all$locs[,1],simulations_all$locs[,2])
  #neigh <- dnearneigh(locs, d1=0, d2=sqrt(2)+10, longlat = FALSE)
  #winnb <- nb2WB(neigh)
  #winnb <- simulations_all$winnb
  #L <- length(winnb$adj)
  #str(winnb)
  #Setting parameters
  data <- list(C,Y, temp, agree, disagree, altitudes)
  #dim_data <- dim(data[[1]])
  #n.cov=ncol(cov)
  n.species <- length(unique(c(data[[1]])))
  n.reported <- length(unique(c(data[[2]])))
  n.sites <- length(data[[1]])
  #alpha <- rep(1, (n.species+1))
  
  # Constants for the MCMC model
  const <- list(n.sites=n.sites,
                n.species=n.species,
                #n.sites = dim_data[3], 
                #n.species = (dim_data[1]), 
                n.reported=n.reported)#, 
                # n.visit=1, #use only for data
                #n.cov=n.cov,
               # L=L, 
               # adj=winnb$adj, 
                #num=winnb$num, 
                #weights=winnb$weights)
  alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.reported))
  for(spe.tag in 1:const$n.species){
    for(k in 1:(const$n.reported)){
      alpha[spe.tag,k] <- 1
    }
  }
  #alpha <- rep(1, (const$n.species+1))
  
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
  
  
  #mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0","beta1","beta2","beta3","beta4","tau","omega", "prop","lambda"))
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
  #GR.diag <- gelman.diag(mcmc.out$samples, multivariate = FALSE)
  #mcmcplot(mcmc.out$samples['beta1',])
  #all(GR.diag$psrf[1:12,"Point est."] < 1.1)
  #which(GR.diag$psrf[1:12,"Point est."] > 1.1) 
all_rhat <- all(Rhat_data < 1.05) #Rhat is less than 1.05
N_over_rhat <-length(which(Rhat_data > 1.05))/length(subset_parameters) #Rhats over 1.04
  
  output <- mcmc.out$summary$all.chains
  beta0 <- output[(1:const$n.species),1]
  beta1 <- output[((const$n.species+1):(const$n.species*2)),1]
  beta2 <- output[(((const$n.species*2)+1):(const$n.species*3)),1]
  beta3 <- output[(((const$n.species*3)+1):(const$n.species*4)),1]
  lambda<- matrix(output[((const$n.species*4)+1) : (((const$n.species*4))+ (const$n.species*(const$n.sites))),1], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  omega <- matrix(output[((((const$n.species*4))+ (const$n.species*(const$n.sites)))+1): ((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))),1], nrow=const$n.species, ncol=(const$n.reported), byrow=FALSE)
  #prop <- matrix(output[(((((const$n.species*2))+ (const$n.species*(const$n.sites)))+(const$n.species*(n.species+1)))) :((((const$n.species*2))+ (const$n.species*(const$n.sites)))+(const$n.species*(n.species+1)))+ (const$n.species*(const$n.sites)),1], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  prop <-matrix(output[((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))+1):((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))+(const$n.species*const$n.sites)),1],nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  tau <- tail(output[,1])[6]
  return(list(omega=omega, beta0=beta0, beta1 = beta1,beta2=beta2,beta3=beta3, prop=prop, tau=tau,lambda=lambda, mcmc.out=mcmc.out,all_rhat = all_rhat,N_over_rhat = N_over_rhat))
}

cl <- makeCluster(6)
setDefaultCluster(cl)
rep_estimates <- pblapply(data_for_nimble,run_simulations, cl=cl)
#rep_estimates <- run_simulations(data_for_nimble)

save(rep_estimates, file="data_results.RData")

