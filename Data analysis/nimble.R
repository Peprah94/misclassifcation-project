#setwd("/Volumes/kwakupa/misclassification/new")
load("nimble_data.RData")
#source("distance_to_roads.R")
library(spatstat)
library(dplyr)
library(sp)
library(gstat)
library(ggplot2)
library(INLA)
library(rgeos)
library(pbapply)
library(nimble)
library(MCMCglmm)
library(coda)
library(MCMCpack)
library(boot)
library(MASS)
library(parallel)
require(ggmcmc)

run_simulations <- function(simulations_all, type, model_selection){
  # packages copied within the code for the parallelisation
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(mcmcplots)
  library(MCMCpack)
  library(boot)
  #library(extraDistr)
  library(MASS)
  #library(spdep)
  require(ggmcmc)
  
  my_proportions <- function(x){
    ret <- proportions(x)
    return(ret)
  }
  
  nimble_proportions <- nimbleRcall(
    prototype = function(
    x=double(1)
    # beta is a vector
    ) {},
    returnType = double(1), # outcome is a vector
    Rfun = 'my_proportions'
  )
  
  start_time <- Sys.time()
  
  
  
  code <- nimbleCode({
    # Priors for beta
    for(spe.tag in 1:(n.species-1)){
      beta0[spe.tag] ~ dnorm(0,sd=10)
      beta1[spe.tag] ~ dnorm(0, sd=10)
      beta2[spe.tag] ~ dnorm(0, sd=10)
     beta3[spe.tag] ~ dnorm(0, sd=10)
     beta4[spe.tag] ~ dnorm(0, sd=10)
     beta5[spe.tag] ~ dnorm(0, sd=10)
    }
    
    
    
    #Classification process prior
    if(type == "variable" | type == "constant" | type == "intercept"|type == "main"|type == "only_principal_cov"){ 
    for(spe.tag in 1:n.species){
      for(report.tag in 1:(n.reported-1)){
        gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 10)
        gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 10)
      }
    }
    }
    
    #user effect
    if( type == "constant" | type == "intercept"|type == "main"| type == "only_principal_cov"|type == "op_cov_finter" ){
    for(user.tag in 1:n.users){
      for(spe.tag in 1:n.species){
      user_effect[spe.tag, user.tag] ~ dnorm(0, sd= 10)
      }
    }
    }
      
    if(type == "variable"){
      for(user.tag in 1:n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:(n.reported-1)){
        user_effect[spe.tag,report.tag, user.tag] ~ dnorm(0, sd= 10)
        }
        }
      }
    }
    

    #main model covariate selection
    if(model_selection == TRUE){
      p ~ dunif(0,1)
      psi ~ dbern(p)
    }else{
      p <- 1
      psi <- 1
}
    
    
    #Use species 1 as reference
    if(type == "variable" | type == "constant" | type == "intercept"| type == "op_cov_finter"| type == "only_principal_cov"){
      for(site.tag in 1:n.sites){
        #for(site.tag in 1:n.sites){
          lambda[1,site.tag] <- 1
          #for(spe.tag in 2:n.species){
          lambda[2, site.tag] <- exp(beta0[1] + beta1[1]*altitude[site.tag] + beta2[1]*distance[site.tag] + beta3[1]*precipitation[site.tag] + beta4[1]*mean_temperature[site.tag])
          lambda[3, site.tag] <- exp(beta0[2] + beta1[2]*altitude[site.tag] + beta2[2]*distance[site.tag] + beta3[2]*precipitation[site.tag] + beta4[2]*mean_temperature[site.tag])
          #}
        #}
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
    }
    
 
    if(type == "variable"){
      
      #Prior for the Confusion Matrix
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
         # for(user.tag in 1:n.users){
          alpha[spe.tag, n.reported, user.tag] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, user.tag] <- exp(gamma0[spe.tag, report.tag] + psi * user_effect[spe.tag, report.tag, user.tag])
           # }
            }
        }
      }
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            #for(user.tag in 1:n.users){
            z.omega[spe.tag, report.tag, user.tag] <- alpha[spe.tag, report.tag, user.tag]/sum(alpha[spe.tag, 1:n.reported, user.tag])
          #  }
}
        }
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag] ])
        }
      }
    }
    
    
    if(type == "intercept"){
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, user.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, user.tag ] <- exp(gamma0[spe.tag, report.tag])
          }
        }
      }
      
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            z.omega[spe.tag, report.tag, user.tag ] <- alpha[spe.tag, report.tag, user.tag ]/sum(alpha[spe.tag, 1:n.reported, user.tag ])
          }
        }
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag]])
        }
      }
    }
    
    
    if(type == "constant"){
      
      for(spe.tag in 1:n.species){
        for(report.tag in 1:n.reported){
          alpha[spe.tag, report.tag] ~ dexp(1)
        }
      }
      # Confusion Matrix for the Species
      for(spe.tag in 1:n.species){
        z.omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }
      
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(z.omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }
    
    if(type == "main"){
      for(site.tag in 1:n.sites){
        #for(site.tag in 1:n.sites){
        lambda[1,site.tag] <- 1
        #for(spe.tag in 2:n.species){
        lambda[2, site.tag] <- exp(beta0[1] + beta1[1]*altitude[site.tag] + beta2[1]*distance[site.tag] + beta3[1]*precipitation[site.tag] + beta4[1]*mean_temperature[site.tag] + psi * user_effect[users[site.tag]])
        lambda[3, site.tag] <- exp(beta0[2] + beta1[2]*altitude[site.tag] + beta2[2]*distance[site.tag] + beta3[2]*precipitation[site.tag] + beta4[2]*mean_temperature[site.tag]+ psi * user_effect[users[site.tag]])
        #}
        #}
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
      
      for(spe.tag in 1:n.species){
        for(report.tag in 1:n.reported){
          alpha[spe.tag, report.tag] ~ dexp(1)
        }
      }
      # Confusion Matrix for the Species
      for(spe.tag in 1:n.species){
        z.omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }
      
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(z.omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }
    
    if(type == "only_principal_cov"){
      
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          # for(user.tag in 1:n.users){
          alpha[spe.tag, n.reported, user.tag] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, user.tag] <- exp(gamma0[spe.tag, report.tag] + psi * user_effect[spe.tag, user.tag])
            # }
          }
        }
      }
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            #for(user.tag in 1:n.users){
            z.omega[spe.tag, report.tag, user.tag] <- alpha[spe.tag, report.tag, user.tag]/sum(alpha[spe.tag, 1:n.reported, user.tag])
            #  }
          }
        }
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag] ])
        }
      }
    }
    
    
    
   if(type == "op_cov_finter"){
      for(spe.tag in 1:n.species){
        gamma0[spe.tag, spe.tag] <- beta5[1]
      }
      gamma0[1, 2]~ dnorm(0, sd = 1)
      gamma0[1, 3]~ dnorm(0, sd = 1)
      gamma0[1, 4]~ dnorm(0, sd = 1)
      gamma0[1, 5]~ dnorm(0, sd = 1)
      
      gamma0[2, 1]~ dnorm(0, sd = 1)
      gamma0[2, 3]~ dnorm(0, sd = 1)
      gamma0[2, 4]~ dnorm(0, sd = 1)
      gamma0[2, 5]~ dnorm(0, sd = 1)
      
      gamma0[3, 1]~ dnorm(0, sd = 1)
      gamma0[3, 2]~ dnorm(0, sd = 1)
      gamma0[3, 4]~ dnorm(0, sd = 1)
      gamma0[3, 5]~ dnorm(0, sd = 1)
      
      #gamma0[1, 2]~ dnorm(0, sd = 10)
      #gamma0[1, 2]~ dnorm(0, sd = 10)
      
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          #gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      #Prior for the Confusion Matrix
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, user.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, user.tag] <- exp(gamma0[spe.tag, report.tag] + psi * user_effect[spe.tag,user.tag])
          }
        }
      }
      
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            z.omega[spe.tag, report.tag, user.tag ] <- alpha[spe.tag, report.tag, user.tag ]/sum(alpha[spe.tag, 1:n.reported, user.tag ])
          }
        }
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag]])
        }
      }
      
    } 
    

    for(validation_site.tag in 1:n.validated){
      ind_all_predictions[validation_site.tag] <- equals(C[1, (index_for_splitting[validation_site.tag])], Validation_C[ validation_site.tag])
    }
    
    for(mismatch_site.tag in 1:n.mismatch){
      ind_all_precisions[mismatch_site.tag] <- equals(C[1, (validation_indices_for_mismatch[mismatch_site.tag])], validation_mismatchC[mismatch_site.tag])
    }
    
    for(correctmatch_site.tag in 1:n.correctmatch){
      ind_all_recall[correctmatch_site.tag] <- equals(C[1, (validation_correct_match[correctmatch_site.tag])], validation_correctmatchC[correctmatch_site.tag])
    }
    
    
    accuracy <- sum(ind_all_predictions[1:n.validated])/ n.validated
    precision <- sum(ind_all_precisions[1:n.mismatch]) / n.mismatch
    recall <- sum(ind_all_recall[1:n.correctmatch])/n.correctmatch
  })
  
  # Retrieving data from the simulated data
  C <- matrix(as.numeric(as.factor(simulations_all$C)), nrow=1) # Unverified citizen science observation
  Y <- matrix(as.numeric(as.factor(simulations_all$Y)), nrow=1) #Verified Citizen science observation
  altitude <- simulations_all$altitude  #Covariates
  distance  <- simulations_all$distance
  users <- simulations_all$users
  Validation_C <- as.numeric(as.factor(simulations_all$Validation_C))
  Validation_Y <- simulations_all$Validation_Y
  validation_indices_for_mismatch <- simulations_all$validation_indices_for_mismatch
  index_for_splitting <- simulations_all$index_for_splitting
  validation_mismatchC <- simulations_all$validation_mismatchC
  precipitation <- simulations_all$precipitation
  mean_temperature <- simulations_all$mean_temperature
  experience <- simulations_all$experience 
  validation_correct_match <- simulations_all$validation_correct_match
  validation_correctmatchC <- simulations_all$validation_correctmatchC
  validation_mismatchC <- as.numeric(sapply(validation_mismatchC, function(x){
    if(x == "Danaus plexippus"){
      ret <- 2
    }
    if(x == "Danaus gilippus"){
    ret <- 1
    }
    if(x == "Limenitis archippus"){
      ret <- 3
    }
    return(ret)
  }))
  
  validation_correctmatchC <- as.numeric(sapply(validation_correctmatchC, function(x){
    if(x == "Danaus plexippus"){
      ret <- 2
    }
    if(x == "Danaus gilippus"){
      ret <- 1
    }
    if(x == "Limenitis archippus"){
      ret <- 3
    }
    return(ret)
  }))
  
  #data <- list(C,Y, altitude,distance, users)
  
  data <- list(C,Y, altitude,distance, Validation_C,Validation_Y,
               validation_indices_for_mismatch, index_for_splitting,
               validation_mismatchC, users, precipitation,
               mean_temperature, experience,
               validation_correctmatchC, validation_correct_match)
  
  #Constants for the model
  dim_data <- dim(data[[1]])
  n.cov=ncol(cov)
  const <- list(n.sites=dim(data[[1]])[2],
                n.species=max(c(data[[1]]), na.rm = TRUE),
                n.reported = max(c(data[[2]]), na.rm = TRUE),
                n.visit=(dim_data[1]), 
                n.cov=n.cov,
                #nref = max(c(data[[1]]))-1,
                n.mismatch = length(data[[7]]),
                n.validated = length(data[[8]]),
                validation_indices_for_mismatch = data[[7]],
                index_for_splitting = data[[8]],
                validation_correct_match = data[[15]],
                n.users = length(unique(data[[13]])),
                                 users = data[[13]],
                n.correctmatch = length(unique(data[[15]]))
                
  )
  
  idm_data <- list(C = data[[1]], 
                   Y = data[[2]], 
                   altitude =data[[3]],
                   distance= data[[4]],
                   Validation_C= data[[5]],
                   Validation_Y = data[[6]],
                   precipitation = data[[11]],
                   mean_temperature = data[[12]],
                   validation_mismatchC = data[[9]],
                   #experience = data[[13]],
                   validation_correctmatchC = data[[14]])
  
  # Initials for the model
  if(type == "variable" | type == "intercept" | type == "op_cov_finter"|type == "only_principal_cov"){
    alpha <- z.omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.users))
    
    for(user.tag in 1:const$n.users){
      for(spe.tag in 1:const$n.species){
        for(report.tag in 1:(const$n.reported)){
          alpha[spe.tag,report.tag, user.tag] <- 1
        }
      }
    }
    
    for(user.tag in 1:(const$n.users)){
      for(spe.tag in 1:(const$n.species)){
        z.omega[spe.tag,, user.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),user.tag])
      }
    }
    
    if(type == "variable"){
      user_effect = array(rnorm(const$n.species*const$n.users*(const$n.reported-1)), dim = c(const$n.species, (const$n.reported-1), const$n.users))
    }else{
      user_effect = matrix(rnorm(const$n.species*const$n.users), nrow = const$n.species, ncol = const$n.users)
    }
    
    idm_inits <- function(){list(beta1 =rep(1, (const$n.species-1)),
                                 z.omega=z.omega,
                                 beta0 = rep(1, (const$n.species-1)),
                                 beta2 = rep(1, (const$n.species-1)),
                                 beta3 = rep(1, (const$n.species-1)),
                                 beta4 = rep(1, (const$n.species-1)),
                                 beta5 = rep(1, (const$n.species-1)),
                                 user_effect = user_effect,
                                 gamma0 = matrix(rep(1,const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = matrix(rep(1,const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported -1)),
                                 alpha = alpha
    )
    }
  }else{
    alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.reported))
    for(spe.tag in 1:const$n.species){
      for(k in 1:(const$n.reported)){
        alpha[spe.tag,k] <- 1
      }
    }
    #alpha <- rep(1, (const$n.species+1))
    
    
    # Initials for the model
    z.omega <- matrix(NA, nrow=const$n.species,
                    ncol = (const$n.reported)) 
    for(spe.tag in 1:(const$n.species)){
      z.omega[spe.tag,]<- rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.reported)])
    }
    
    idm_inits <- function(){list(beta1 =rep(1, (const$n.species-1)),
                                 z.omega=z.omega,
                                 beta2 = rep(1, (const$n.species-1)),
                                 beta3 = rep(1, (const$n.species-1)),
                                 beta4 = rep(1, (const$n.species-1)),
                                 beta0 = rep(1, (const$n.species-1)),
                                 beta5 = rep(1, (const$n.species-1)),
                                 user_effect = rnorm(const$n.users),
                                 gamma0 = matrix(rep(1,const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = matrix(rep(1,const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported -1)),
                                 alpha = alpha
    )
    }
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
  #nimbleOptions(stop_after_processing_model_code = TRUE)
  mwtc <- nimbleModel(code,data = idm_data, 
                      constants = const, 
                      inits = initsList,
                      dimensions =list(lambda = c(const$n.species, const$n.sites), 
                                       calculate = FALSE))
  
  # Create the model in C
  Cmwtc <- compileNimble(mwtc, 
                         showCompilerOutput = FALSE)
  
  
  
  mcmcconf <- configureMCMC(Cmwtc, 
                            monitors = c("beta0","beta1", "beta2", "beta3", "beta4","beta5",
                                         "gamma0", "gamma1",
                                         "accuracy", "precision", "p", 
                                         "user_effect", "recall", "z.omega"))
  
  #mcmcconf$removeSamplers(c("omega[1,1:4]","omega[2,1:4]", "omega[3,1:4]"))
  #mcmcconf$addSampler(c("omega[1,1:2]"), "RW_dirichlet", adaptive=TRUE, scale=3)
  #mcmcconf$addSampler(c("omega[2,1:2]"), "RW_dirichlet", adaptive=TRUE, scale=3)
  
  Rmcmc <- buildMCMC(mcmcconf)
  
  # Compile 
  cmcmc <- compileNimble(Rmcmc, 
                         project = Cmwtc,
                         resetFunctions = TRUE)
  
  #ret <- cmcmc$run(niter = 1000)
  # Run the MCMC
  mcmc.out <- runMCMC(cmcmc, 
                      niter = 10000,
                      nchains =3,
                      nburnin = 5000,
                      inits = initsList,
                      thin = 5, 
                      setSeed = TRUE, 
                      samples=TRUE, 
                      samplesAsCodaMCMC = TRUE, 
                      summary = TRUE, 
                      WAIC = FALSE)
  
  # Return the summary of MCMC results
  output <- mcmc.out$summary$all.chains
  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  
  mcmclist <- ggs(mcmc.out$samples)
  aa <- mcmclist%>%
    filter(grepl('gamma|beta|user_effect',Parameter))
  
  subset_parameters<- unique(aa$Parameter)
  
  # Rhat values
  subset_Rhat <- aa%>%
    ggs_Rhat()
  Rhat_data <- subset_Rhat$data[,5]
  all_rhat <- all(Rhat_data < 1.1) #Rhat is less than 1.05
  N_over_rhat <-length(which(Rhat_data > 1.1))/length(subset_parameters) #Rhats over 1.04
  
  accuracy <- output[1,2]
  beta0 <- output[2:3,1]
  beta1 <- output[4:5,1]
  beta2 <- output[6:7,1]
  beta3 <- output[8:9,1]
  beta4 <- output[10:11,1]
  beta5 <- output[12:13,1]
  #lambda<- matrix(output[((const$n.species*2)+1) : (((const$n.species*2))+ (const$n.species*(const$n.sites))),2], nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  gamma0 <- matrix(output[14:28,1], nrow=const$n.species, ncol=(const$n.reported-1), byrow=FALSE)
  #gamma0 <- matrix(output[4:5,1], nrow=(const$n.species), ncol=(const$n.reported-1), byrow=FALSE)
  #gamma1 <- matrix(output[8:11,1], nrow= 2, ncol= 2, byrow=FALSE)
  gamma1 <- matrix(output[29:43,1], nrow=(const$n.species), ncol=(const$n.reported-1), byrow=FALSE)
  p <- output[44,1]
  precision <- output[45,2]
  recall <- output[46,2]
  if(type != "variable"){
user_effect <- matrix(output[47:64,1], nrow = const$n.species, ncol = const$n.users, byrow = FALSE)
  }else{
    user_effect <- array(output[47:136,1], dim = c(const$n.species, (const$n.reported-1), const$n.users))
  }
  #omega <- array(output[53:160,1], dim = c(const$n.sites, const$n.reported, const$n.users), byrow = FALSE)
  #user_effect <- output[26:37,1]
  #contrasts <- matrix(tail(output[,2], const$n.species*const$n.species), nrow=const$n.species, ncol=(const$n.species), byrow=FALSE)
  #prop <-matrix(tail(output[,2],(const$n.species*(const$n.sites))),nrow=const$n.species, ncol=const$n.sites, byrow=FALSE)
  return(list(beta0=beta0, beta1 = beta1,beta2 = beta2,beta3 = beta3,
              beta4 = beta4, beta5 = beta5,gamma0= gamma0, gamma1 = gamma1,  
              accuracy= accuracy, precision = precision, p = p, 
              user_effect = user_effect, recall = recall,output = output 
              ))
}



# rep_estimates <- lapply(list( "main","variable", "intercept", "constant"), function(x){
#   ret <- run_simulations(dataout, x)
# }
# )
# 
# save(rep_estimates, file="estimated_data.RData")
