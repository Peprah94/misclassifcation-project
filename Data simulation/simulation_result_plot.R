# Plots results from the simulation study 
# Takes inputs from the data returned by "estimate_sim.R" and the simulated data

#Packages needed to run the script
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggridges)
library(grid)
library(spatstat) 

#Load data
load("simulated data1.RData")
load("estimateddata_new.RData")

#True parameters used for simulation
input <- list(
  constants=list(
    n.species=3,
    #n.sites=nrow(locs),
    n.visit=5, 
    p.tag= 0.7,
    dim = c(20,20),
    sample.size = 400
  ),
  plotdata = list(
    log_lambda=FALSE,
    prop_plot=FALSE),
  fixed_effects = list(
    beta0 = c(-1,1,0),
    beta1=c(2,0,-2),
    beta2=c(0.8,-1,0.5)
  ),
  hyperparameters = list(
    sigma2x = c(0.5, 0.5, 0.5),
    kappa = c(1.5,1.5, 1.5)
  ),
  covariates = list(
    cov = cov
  )
)


cov_test <- seq(-1,1,0.01)
cov_test <- rbind(cov_test,cov_test)
n.species=3
lambda_true <-prop_true <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
max_lambda_true <- max_prob_true <- vector("numeric", length(cov_test[1,]))

for(site.tag in 1:length(cov_test[1,])){
  for(spe.tag in 1:n.species){
    #True and estimated abundance
    lambda_true[spe.tag, site.tag] <- exp(input$fixed_effects$beta0[spe.tag] + input$fixed_effects$beta1[spe.tag]*cov_test[1,site.tag])
    #estimation proportions
    prop_true[, site.tag] <- proportions(lambda_true[1:n.species, site.tag])
  }
  max_lambda_true[site.tag] <- which.max(lambda_true[,site.tag])
}

#Put data together
data_all <- data.frame(cov_test=cov_test[1,], 
                       true_lambda= as.factor(max_lambda_true),
                       test_species1 = lambda_true[1,],
                       test_species2 = lambda_true[2,],
                       test_species3 = lambda_true[3,],
                       true_state1 = prop_true[1,],
                       true_state2 = prop_true[2,],
                       true_state3 = prop_true[3,])

#Plot the covariate range and the predicted probability
# FIgure 1 in main paper
data_all%>%
  dplyr:: select(cov_test, true_state1,true_state2,true_state3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot(aes(x=cov_test, y=value,  fill= variable))+
  geom_density_line(stat = "identity", size=0.5, alpha=0.5, position = "fill")+
  scale_fill_manual(name='', values=c("true_state1" = "green4", "true_state2" = "red", "true_state3" = "blue"))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("Covariate range")+
  ylab(expression(Probability (p[i][j])))+
  theme(legend.title = element_blank(), axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
ggsave("covariates.png")

###############
# Estimates of Parameters
################
#Function to estimate RMSE
rmse <- function(x){
  sqrt(mean(x^2))
}
###################
#contrast of betas
#############
contrast_beta <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    x[[4]]
  }
})
contrast_beta_dataframe <- do.call("rbind", contrast_intercept)
median_beta <- apply(contrast_beta_dataframe,2,median)

#coverage of contasts of betas
coverage <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    x[[5]]
  }
})
coverage_dataframe <- do.call("rbind", coverage) 
mean_coverage <- apply(coverage_dataframe,2,mean)

# bias in contrasts
bias_contrasts <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    c(x[[4]] - c(2,1,-2,-4))
  }
})
bias_contrasts_dataframe <- do.call("rbind", bias_contrasts)
#rmse of contrasts
rmse_contrasts <- apply(bias_contrasts_dataframe,2,rmse)

##################################
#classification probabilities (omega)
##################################

bias_confusion <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    c(x[[1]] - sim[[1]]$omega)
  }
})
bias_confusion_dataframe <- do.call("rbind", bias_confusion)
#rmse of classification probabilities
rmse_confusion <- apply(bias_confusion_dataframe,2,rmse)

###########################
#Testing the results on newdata
#############################
#simulating new covariate
cov_test <- runif(600,-1,1)
cov_test <- rbind(cov_test, cov_test)

#input values for simulating new data
input_test <- list(
  constants=list(
    n.species=3,
    #n.sites=nrow(locs),
    n.visit=5, 
    p.tag= c(0.7,0.8,0.9),
    dim = c(30,20),
    sample.size = 600
  ),
  plotdata = list(
    log_lambda=FALSE,
    prop_plot=FALSE),
  fixed_effects = list(
    beta0 = c(-1,1,0),
    beta1=c(2,0,-2),
    beta2=c(0.8,-1,0.5)
  ),
  hyperparameters = list(
    sigma2x = c(0.5, 0.5, 0.5),
    kappa = c(1.5,1.5, 1.5)
  ),
  covariates = list(
    cov = cov_test
  )
)

dim = c(30,20)
win <- owin(c(0,dim[1]), c(0,dim[2])) # keeping this rectangular is really important to check for errors in the code - otherwise easy to get x and y confused
# set number of pixels
spatstat.options(npixel=c(dim[1],dim[2]))
y0 <- seq(win$yrange[1], win$yrange[2],
          length=spatstat.options()$npixel[2])
x0 <- seq(win$xrange[1], win$xrange[2],
          length=spatstat.options()$npixel[1])
locs <- expand.grid(x0,y0)
sim_test <- genData(input_test, seed = 1000) # From the spatial_sim.R function

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

#Confusion matrix
omega <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    c(x[[1]])
  }
})
combine_omega <- do.call("rbind", omega)

#median of confusion matrix
median_omega <- apply(combine_omega, 2, median)
omega_est <- matrix(median_omega, nrow = 3, ncol=4, byrow=F)

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

library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)


g3 <- test_results(1)%>%
  dplyr:: select(cov_test, species1,species2,species3)%>%
  #dplyr::mutate(new_species2 = species1+ species2,
  #              new_species3=species1+ species2+species3)%>%
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

# Plots of results in figure 5
gg_sub <- ggarrange(g3,g6,g9,g12, ncol=2, nrow=2, common.legend = TRUE, legend="right", labels = c("a)", "b)", "c)", "d)"))
ggpubr::annotate_figure(gg_sub, left = text_grob("Verified species probability from estimated parameters", rot = 90, vjust = 1),
                        bottom = text_grob("Covariate values"))
ggsave("results_sub.png")
