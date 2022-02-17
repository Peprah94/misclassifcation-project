# packages needed to run this script
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggridges)
library(grid)

source("simulations.R")

#Load data from simulation.R script
load("simulated data1.RData")

#Load data from estimate_sim.R script
load("estimateddata_new.RData")

#True parameters used in the simuation (same as used in simulation.R script)
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

######################################
# Code for figure 1 in paper
################################
cov_test <- seq(-2,2,0.01) #covariate range to plot
n.species=3 #number of species

#Estimating the abundance and proportion for this covariate value
lambda_true <-prop_true <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
  for(site.tag in 1:length(cov_test[1,])){
    for(spe.tag in 1:n.species){
      #True and estimated abundance
      lambda_true[spe.tag, site.tag] <- exp(input$fixed_effects$beta0[spe.tag] + input$fixed_effects$beta1[spe.tag]*cov_test[site.tag])

      #estimation proportions
      prop_true[, site.tag] <- proportions(lambda_true[1:n.species, site.tag])

    }
    max_lambda_true[site.tag] <- which.max(lambda_true[,site.tag])
  }

# Return dataframe
data_all <- data.frame(cov_test=cov_test[1,], 
                   true_state1 = prop_true[1,],
                   true_state2 = prop_true[2,],
                   true_state3 = prop_true[3,])

# Plotting the data
data_all%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot(aes(x=cov_test, y=value,  fill= variable))+
  #geom_line(aes(y=cov_test, x=value, col=variable))+
  geom_density_line(stat = "identity", size=0.5, alpha=0.5, position = "fill")+
  scale_fill_manual(name='', values=c("true_state1" = "green4", "true_state2" = "red", "true_state3" = "blue"))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("")+
  ylab(expression(Probability (p[i][j])))+
  theme(legend.title = element_blank())

#saving the plot
ggsave("covariates.png")

##################################
# Extracting parameters from the results
#################################
#Extracting the median of all the estimates of beta0 and beta
beta0 <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    x[[2]] 
  }
})
beta0_dataframe <- do.call("rbind", beta0)

beta1 <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    x[[3]] 
  }
})
beta1_dataframe <- do.call("rbind", beta1)

#estimating median of beta0 and beta1
combine_beta <- cbind(beta0_dataframe, beta1_dataframe)
median_beta <- apply(combine_beta,2,median)
median_estimates <- list(
  beta0 = median_beta[1:3],
  beta1=median_beta[4:6]
)
quantiles_beta <- apply(combine_beta, 2, function(x){quantile(x, c(0.05,0.975))})

#estimating the precision of the estimates
sd_beta <- apply(combine_beta,2,sd)
precision_beta <- 1/sd_beta

#Confusion matrix
omega <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    c(x[[1]])
  }
})
combine_omega <- do.call("rbind", omega)

#median and precision of confusion matrix
median_omega <- apply(combine_omega, 2, median)
omega_est <- matrix(median_omega, nrow = 3, ncol=4, byrow=F)
quantiles_omega <- apply(combine_omega, 2, function(x){quantile(x, c(0.05,0.975))})

sd_omega <- apply(combine_omega, 2, sd)
precision_omega <- 1/sd_omega

print(c(precision_beta, precision_omega))

###############################################
# Plotting data for results from simulations
##############################################

#estimates for the bias of constrast in intercept parameters (baseline reference is species 1)
contrast_intercept <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    z <- outer(x[[2]],x[[2]],'-')
    new_z <- z[lower.tri(z)]
    
    z1 = outer(input$fixed_effects$beta0,input$fixed_effects$beta0,'-') 
    new_z1 <- z1[lower.tri(z1)]
    
    (new_z - new_z1)[1:2]
    }
})
contrast_intercept_dataframe <- do.call("rbind", contrast_intercept)


#estimates for the bias of constrast in covariate parameters (baseline reference is species 1)
bias_covariate <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    z <- outer(x[[3]],x[[3]],'-')
    new_z <- z[lower.tri(z)]
    
    z1 = outer(input$fixed_effects$beta1,input$fixed_effects$beta1,'-') 
    new_z1 <- z1[lower.tri(z1)]
    
    (new_z - new_z1)[1:2]
#x[[3]] - input$fixed_effects$beta1
  }
})
bias_covariate_dataframe <- do.call("rbind", bias_covariate)

#estimates for the bias of classification/confusion matrix
bias_confusion <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    c(x[[1]] - sim[[1]]$omega)
  }
})
bias_confusion_dataframe <- do.call("rbind", bias_confusion)

#Plotting the contrasts

all_data <- data.frame(cbind(contrast_intercept_dataframe,bias_covariate_dataframe))
colnames(all_data) <- c("beta21","beta31",
                       "beta11", "beta12")

all_data_new <- melt(all_data)

ggplot(all_data_new)+
  geom_boxplot(aes(value, variable))+
  theme(legend.title = "none")+
  geom_vline(xintercept=0, linetype="dotted",color="red", size=1)+
   scale_y_discrete(labels=c('beta21' = expression(Delta*(beta[0][2] - beta[0][1])),
                           'beta31' = expression(Delta*(beta[0][3] - beta[0][1])),
                            'beta11'=expression(Delta*(beta[1][2] - beta[1][1])),
                      'beta12'=expression(Delta*(beta[1][3] - beta[1][1]))
  ))+
  theme_bw()+
  xlab("Bias in estimates")+
  ylab("Parameters")

#save the plot
ggsave("bias_boxplot.png")

#############################
# Predicting true values for the new data
##############################

#New covariate simulated
cov_test <- seq(-3,2.99,0.01)

input_test <- list(
  constants=list(
    n.species=3,
    n.sites=600,
    n.visit=5, 
    p.tag= c(0.7,0.8,0.9)
  ),
  fixed_effects = list(
    beta0 = c(-3,1,2),
    beta1=c(1,0,-1),
  ),
  covariates = list(
    cov = cov_test
  )
)

# Simulating new data 
sim_test <- genData(input_test, seed = 1000) # From the simulations.R function

#prediction of true species
n.species=3
lambda_true <- lambda_test <- prop_true <- prop_test <-prob_return_true <- prob_return_false <- prob_ret_true <- prob_ret_false <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
max_lambda_true <- max_lambda_test <- max_prob_true <- max_prob_test <- vector("numeric", length(cov_test[1,]))
omega_true = sim_test$omega

test_results <- function(category){
  lambda_true <- lambda_test <- prop_true <- prop_test <-prob_return_true <- prob_return_false <- prob_ret_true <- prob_ret_false <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
  max_lambda_true <- max_lambda_test <- max_prob_true <- max_prob_test <- vector("numeric", length(cov_test[1,]))
reported_species = rep(category, length(cov_test[1,]))
for(site.tag in 1:length(cov_test[1,])){
  for(spe.tag in 1:n.species){
    #True and estimated abundance
    lambda_true[spe.tag, site.tag] <- exp(input$fixed_effects$beta0[spe.tag] + input$fixed_effects$beta1[spe.tag]*cov_test[1,site.tag])
    lambda_test[spe.tag,site.tag] <- exp(median_estimates$beta0[spe.tag]+ median_estimates$beta1[spe.tag]*cov_test[1,site.tag])
    #estimation proportions
    prop_true[, site.tag] <- proportions(lambda_true[1:n.species, site.tag])
    prop_test[,site.tag] <- proportions(lambda_test[1:n.species, site.tag])
  }
 }
    for(site.tag in 1:length(cov_test[1,])){
     for(spe.tag in 1:n.species){
    #Estimating return probs
    prob_return_true[spe.tag, site.tag] <- prop_true[spe.tag, site.tag]*omega_true[spe.tag,reported_species[site.tag]]
    prob_return_false[spe.tag, site.tag] <- prop_test[spe.tag, site.tag]*omega_est[spe.tag,reported_species[site.tag]]
     }
    }
    for(site.tag in 1:length(cov_test[1,])){
      for(spe.tag in 1:n.species){
    # Sum over the species
    prob_ret_true[spe.tag, site.tag] <- prob_return_true[spe.tag, site.tag]/(colSums(prob_return_true)[site.tag])
    prob_ret_false[spe.tag, site.tag] <- prob_return_false[spe.tag, site.tag]/(colSums(prob_return_false)[site.tag])
      }
      }

#Return the species with maximum intensity and probability
for(site.tag in 1:length(cov_test[1,])){
  max_lambda_true[site.tag] <- which.max(lambda_true[,site.tag])
  max_lambda_test[site.tag] <- which.max(lambda_test[,site.tag])
  max_prob_true[site.tag] <- which.max(prob_ret_true[,site.tag])
  max_prob_test[site.tag] <- which.max(prob_ret_false[,site.tag])
}

data <- data.frame(cov_test=cov_test[1,], 
                   true_lambda= as.factor(max_lambda_true),
                   test_lambda= as.factor(max_lambda_test),
                   true_prob = as.factor(max_prob_true),
                   test_prob = as.factor(max_prob_test),
                   species1 = prob_ret_false[1,],
                   species2 = prob_ret_false[2,],
                   species3 = prob_ret_false[3,],
                   true_species1 = prop_true[1,],
                   true_species2 = prop_true[2,],
                   true_species3 = prop_true[3,])
return(data)
}

#####
#reported species = species 1
g3 <- test_results(1)%>%
  dplyr:: select(cov_test, species1,species2,species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot(aes(x=cov_test, y=value, fill= variable))+
  #geom_line(aes(y=cov_test, x=value, col=variable))+
  geom_density_line(stat = "identity", size=0.5, alpha=0.5, position = "fill")+
  scale_fill_manual(name='', values=c("species1" = "green4", "species2" = "red", "species3" = "blue"))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank())+
  xlab("")

##########
#reported species = species 2
g6 <- test_results(2)%>%
  dplyr:: select(cov_test, species1,species2,species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot(aes(x=cov_test, y=value,  fill= variable))+
  #geom_line(aes(y=cov_test, x=value, col=variable))+
  geom_density_line(stat = "identity", size=0.5, alpha=0.5, position = "fill")+
  scale_fill_manual(name='', values=c("species1" = "green4", "species2" = "red", "species3" = "blue"))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("")+
  ylab("Verified species prob \n from estimated parameters")+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")

####################
#reported species = species 3
g9 <- test_results(3)%>%
  dplyr:: select(cov_test, species1,species2,species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot(aes(x=cov_test, y=value,  fill= variable))+
  #geom_line(aes(y=cov_test, x=value, col=variable))+
  geom_density_line(stat = "identity", size=0.5, alpha=0.5, position = "fill")+
  scale_fill_manual(name='', values=c("species1" = "green4", "species2" = "red", "species3" = "blue"))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("")+
  ylab("Verified species prob \n from estimated parameters")+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")

#######
#reported species = species 4
g12 <- test_results(4)%>%
  dplyr:: select(cov_test, species1,species2,species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot(aes(x=cov_test, y=value,  fill= variable))+
  #geom_line(aes(y=cov_test, x=value, col=variable))+
  geom_density_line(stat = "identity", size=0.5, alpha=0.5, position = "fill")+
  scale_fill_manual(name='', values=c("species1" = "green4", "species2" = "red", "species3" = "blue"))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("")+
  ylab("Verified species prob \n from estimated parameters")+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")

# Plot of figure 4
gg_sub <- ggarrange(g3,g6,g9,g12, ncol=2, nrow=2, common.legend = TRUE, legend="bottom", labels = c("a)", "b)", "c)", "d)"))
ggpubr::annotate_figure(gg_sub, left = text_grob("Verified species probability from estimated parameters", rot = 90, vjust = 1),
                        bottom = text_grob("Covariate values"))
ggsave("results_sub.png")

