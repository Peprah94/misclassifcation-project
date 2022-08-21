#This code is used to estimate the parameters for the model
#for the simulation scenarios


# Load packages needed to run the code
library(pbapply)
library(nimble)
library(MCMCglmm)
library(coda)
library(MCMCpack)
library(boot)
library(MASS)
library(parallel)
require(ggmcmc)
op <- pboptions(type="timer")

# Function that estimates the parameters
#Needed to run the code in parallel
run_simulations <- function(simulations_all, type, model_selection){
  # packages copied within the code for the parallelisation
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(mcmcplots)
  library(MCMCpack)
  library(boot)
  library(MASS)
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
    }
    
    #main model covariate selection
    if(model_selection == TRUE){
    p ~ dunif(0,1)
    psi ~ dbern(p)
    }
    
    if(model_selection == FALSE){
      p <- 1
      psi <- 1
    }

    #Use species 1 as reference
    if(type == "variable" | type == "constant" | type == "intercept"){
      #Prior for the classification process
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      #latent mean intensity of the ecological state
      for(site.tag in 1:n.sites){
        lambda[1,site.tag] <- 1
        lambda[2, site.tag] <- exp(beta0[1] + beta1[1]*cov[site.tag])
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
    }
    
    
    if(type == "variable"){
 
    #Prior for the Confusion Matrix
    for(site.tag in 1: n.sites){
    for(spe.tag in 1:n.species){
      alpha[spe.tag, n.reported, site.tag ] <- 1
      for(report.tag in 1:(n.reported-1)){
        alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag, report.tag] * cov_omega[site.tag])
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
      
      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      
      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
        }
      }
    }
    
    
    if(type == "intercept"){
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
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
      
      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      
      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
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
        omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }
      
      for(visit.tag in 1:n.visit){
      for(site.tag in 1:n.sites){
        C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
      for(site.tag in 1:n.sites){
        Y[visit.tag, site.tag] ~ dcat(omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }
    
    if(type == "main"){
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      for(site.tag in 1:n.sites){
        lambda[1,site.tag] <- 1
        lambda[2, site.tag] <- exp(beta0[1] + beta1[1]*cov[site.tag]) + beta2[1]*cov_omega[site.tag]
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
        omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }
      
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }

    if(type == "only_principal_cov"){
      #Mean intensity
      for(site.tag in 1:n.sites){
        lambda[1,site.tag] <- 1
        lambda[2, site.tag] <- exp(beta0[1] + beta1[1]*cov[site.tag])
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
    
    #prior for the classification process
      for(spe.tag in 1:n.species){
        gamma1[spe.tag] ~ dnorm(0, sd = 1)
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }

      
      #Prior for the Confusion Matrix
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag] * cov_omega[site.tag])
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
      
      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      
      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
        }
      }
      
    }

    
    if(type == "op_cov_finter"){
      #Mean intensity
      for(site.tag in 1:n.sites){
        lambda[1,site.tag] <- 1
        lambda[2, site.tag] <- exp(beta0[1] + beta1[1]*cov[site.tag])
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
      
      #Prior for the classification process
      for(spe.tag in 1:n.species){
        gamma1[spe.tag] ~ dnorm(0, sd = 1)
          gamma0[spe.tag, spe.tag] <- beta2[1]
      }
      gamma0[1, 2]~ dnorm(0, sd = 1)
      gamma0[2, 1]~ dnorm(0, sd = 1)
      
      #Prior for the Confusion Matrix
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag] * cov_omega[site.tag])
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
      
      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      
      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
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
  C <- simulations_all$C # verified state
  Y <- simulations_all$Y #Reported category
  Validation_C <- simulations_all$Validation_C
  Validation_Y <- simulations_all$Validation_Y
  validation_indices_for_mismatch <- simulations_all$validation_indices_for_mismatch
  index_for_splitting <- simulations_all$index_for_splitting
  validation_mismatchC <- simulations_all$validation_mismatchC
  cov <- simulations_all$cov #Covariates
  cov_omega <- simulations_all$cov_omega
  validation_correct_match <- simulations_all$validation_indices_for_match
  validation_correctmatchC <- simulations_all$validation_matchC
  data <- list(C,Y, cov,cov_omega, Validation_C,Validation_Y,
               validation_indices_for_mismatch, index_for_splitting,
               validation_mismatchC,  validation_correct_match,
               validation_correctmatchC)
  
  #Constants for the model
  dim_data <- dim(data[[1]])
  n.cov=ncol(cov)
  const <- list(n.sites=dim_data[2],
                n.species=max(c(data[[1]]), na.rm = TRUE),
                n.reported = max(c(data[[2]]), na.rm = TRUE),
                n.visit=(dim_data[1]), 
                n.cov=n.cov,
                n.mismatch = length(data[[7]]),
                n.validated = length(data[[8]]),
                n.correctmatch = length(data[[10]]),
                validation_indices_for_mismatch = data[[7]],
                index_for_splitting = data[[8]],
                validation_correct_match = data[[10]]
  )
  
  # Data for the Model
  idm_data <- list(C = data[[1]], 
                   Y = data[[2]], 
                   cov=data[[3]],
                   cov_omega = data[[4]],
                   Validation_C= data[[5]],
                   Validation_Y = data[[6]],
                   validation_mismatchC = data[[9]],
                   validation_correctmatchC = data[[11]])
  
  # Initials for the model
  if(type == "variable" | type == "intercept"){
  alpha <- omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.sites))
  
  for(site.tag in 1:const$n.sites){
  for(spe.tag in 1:const$n.species){
    for(report.tag in 1:(const$n.reported)){
      alpha[spe.tag,report.tag, site.tag] <- 1
    }
  }
  }
    
  for(site.tag in 1:(const$n.sites)){
  for(spe.tag in 1:(const$n.species)){
    omega[spe.tag,, site.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),site.tag])
  }
  }
  
  idm_inits <- function(){list(beta1 =rnorm((const$n.species-1)),
                               omega=omega,
                               beta0 = rnorm((const$n.species-1)),
                               gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                               gamma1 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported -1)),
                               alpha = alpha
                               )
  }
  }else{if(type =="only_principal_cov"| type == "op_cov_finter"){
    alpha <- omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.sites))
    
    for(site.tag in 1:const$n.sites){
      for(spe.tag in 1:const$n.species){
        for(report.tag in 1:(const$n.reported)){
          alpha[spe.tag,report.tag, site.tag] <- 1
        }
      }
    }
    
    for(site.tag in 1:(const$n.sites)){
      for(spe.tag in 1:(const$n.species)){
        omega[spe.tag,, site.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),site.tag])
      }
    }
    
    idm_inits <- function(){list(beta1 =rnorm((const$n.species-1)),
                                 omega=omega,
                                 beta0 = rnorm((const$n.species-1)),
                                 gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = rnorm(const$n.species),
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
    omega <- matrix(NA, nrow=const$n.species,
                    ncol = (const$n.reported)) 
    for(spe.tag in 1:(const$n.species)){
      omega[spe.tag,]<- rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.reported)])
    }
    
    idm_inits <- function(){list(beta1 =rep(1, (const$n.species-1)),
                                 omega=omega,
                                 beta0 = rep(1, (const$n.species-1)),
                                 gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported -1)),
                                 alpha = alpha
    )
    }
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
  mwtc <- nimbleModel(code,data = idm_data, 
                      constants = const, 
                      inits = initsList,
                      dimensions =list(lambda = c(const$n.species, const$n.sites), 
                                       calculate = FALSE))
  
  # Create the model in C
  Cmwtc <- compileNimble(mwtc, 
                         showCompilerOutput = FALSE)
  
  
  
  mcmcconf <- configureMCMC(Cmwtc, 
                            monitors = c("beta0","beta1",
                                         "gamma0","gamma1",
                                         "accuracy", "precision", "p", "recall"))
  
  Rmcmc <- buildMCMC(mcmcconf)
  
  # Compile 
  cmcmc <- compileNimble(Rmcmc, 
                         project = Cmwtc,
                         resetFunctions = TRUE)

  # Run the MCMC
  mcmc.out <- runMCMC(cmcmc, 
                      niter = 10000,
                      nchains =3,
                      nburnin = 5000,
                      inits = initsList,
                      #thin = 5, 
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
    filter(grepl('gamma|beta',Parameter))
  
  subset_parameters<- unique(aa$Parameter)
  
  # Rhat values
  subset_Rhat <- aa%>%
    ggs_Rhat()
  Rhat_data <- subset_Rhat$data[,5]
  all_rhat <- all(Rhat_data < 1.1) #Rhat is less than 1.05
  N_over_rhat <-length(which(Rhat_data > 1.1))/length(subset_parameters) #Rhats over 1.04
  
  #Coverage of contrasts
  coverage <- function(output){
    ifelse(output[4] <= output[6] & output[5] >= output[6],1,0)
  }
  
  true_values <- c( -2, -2,
                   simulations_all$true_values)
  
  df <- cbind(output[2:11,], true_values) #11
  coveragevalues = apply(df,1, coverage)
  
  returnedcoverages <- c(coveragevalues,N_over_rhat=N_over_rhat, all_rhat = all_rhat)
  
  accuracy <- output[1,2]
  beta0 <- output[2,1]
  beta1 <- output[3,1]
  gamma0 <- matrix(output[4:7,1], nrow=2, ncol=(const$n.reported-1), byrow=FALSE)
  if(type == "main" |type == "constant"| type == "intercept" | type == "variable"){
  gamma1 <- matrix(output[8:11,1], nrow= 2, ncol= 2, byrow=FALSE)
  p <- output[12,2]
  precision <- output[13,2]
  recall <- output[14,2]
  }else{
    gamma1 <- matrix(output[8:9,1], nrow= 1, ncol= 2, byrow=FALSE)
    p <- output[10,2]
    precision <- output[11,2] 
    recall <- output[12,2]
  }

  return(list(beta0=beta0, beta1 = beta1, gamma0= gamma0, gamma1 = gamma1, returnedcovarages=returnedcoverages, accuracy= accuracy, precision = precision, p = p, recall = recall))
}

