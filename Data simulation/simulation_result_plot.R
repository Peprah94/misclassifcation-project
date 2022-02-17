library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggridges)

#Load data
load("/Volumes/kwakupa/misclass1/simulated data1.RData")

load("/Volumes/kwakupa/misclass1/estimateddata_new.RData")
rep_est1 <- rep_estimates

load("/Volumes/kwakupa/missclass2/estimateddata_new.RData")
rep_est2 <- rep_estimates

load("/Volumes/kwakupa/missclass3/estimateddata_new.RData")
rep_est3 <- rep_estimates

load("/Volumes/kwakupa/missclass4/estimateddata_new.RData")
rep_est4 <- rep_estimates

load("/Volumes/kwakupa/missclass5/estimateddata_new.RData")
rep_est5 <- rep_estimates

load("/Volumes/kwakupa/missclass6/estimateddata_new.RData")
rep_est6 <- rep_estimates

load("/Volumes/kwakupa/misclass7/estimateddata_new.RData")
rep_est7 <- rep_estimates

load("/Volumes/kwakupa/misclass8/estimateddata_new.RData")
rep_est8 <- rep_estimates

load("/Volumes/kwakupa/misclass9/estimateddata_new.RData")
rep_est9 <- rep_estimates

#load("/Volumes/kwakupa/misclass10/estimateddata_new.RData")
#rep_est10 <- rep_estimates

#True parameters
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
    beta0 = c(-3,1,2),
    beta1=c(1,0,-1),
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


cov_test <- seq(-2,2,0.01)
cov_test <- rbind(cov_test,cov_test)
n.species=3
lambda_true <-prop_true <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
max_lambda_true <- max_prob_true <- vector("numeric", length(cov_test[1,]))
#sim_test <- genData(input)
#omega_true = sim_test$omega
#omega_est = rep_estimates$omega
#reported_species = sim_test$Y[1,]

  #lambda_true <- prop_true <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
  #max_lambda_true <- max_lambda_test <- max_prob_true <- max_prob_test <- vector("numeric", length(cov_test[1,]))
  #reported_species = rep(category, length(cov_test[1,]))
  for(site.tag in 1:length(cov_test[1,])){
    for(spe.tag in 1:n.species){
      #True and estimated abundance
      lambda_true[spe.tag, site.tag] <- exp(input$fixed_effects$beta0[spe.tag] + input$fixed_effects$beta1[spe.tag]*cov_test[1,site.tag])
      #lambda_test[spe.tag,site.tag] <- exp(median_estimates$beta0[spe.tag]+ median_estimates$beta1[spe.tag]*cov_test[1,site.tag])
      #prop_true[spe.tag, site.tag] <- proportions(lambda_true[1:n.species, site.tag])
      #estimation proportions
      prop_true[, site.tag] <- proportions(lambda_true[1:n.species, site.tag])
      #prop_test[,site.tag] <- proportions(lambda_test[1:n.species, site.tag])
    }
    max_lambda_true[site.tag] <- which.max(lambda_true[,site.tag])
  }

data_all <- data.frame(cov_test=cov_test[1,], 
                   true_lambda= as.factor(max_lambda_true),
                   test_species1 = lambda_true[1,],
                   test_species2 = lambda_true[2,],
                   test_species3 = lambda_true[3,],
                   true_state1 = prop_true[1,],
                   true_state2 = prop_true[2,],
                   true_state3 = prop_true[3,])

data_all%>%
  dplyr:: select(cov_test, true_state1,true_state2,true_state3)%>%
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

ggsave("covariates.png")

k3 <- data_all%>%
  dplyr:: select(cov_test, test_species1,test_species2,test_species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot()+
  geom_line(aes(x=cov_test, y=log(value), col=variable))+
  theme_bw()+
  xlab("Covariate")+
  ylab("Log Intensity")+
  theme(legend.title = element_blank())
ggarrange(k1,k3, nrow=1, ncol=2, common.legend = T, legend = "bottom")



# All data
#rep_estimates <- c(rep_est1, rep_est2, rep_est3, rep_est4,rep_est5) #change later
rep_estimates <- c(rep_est1, rep_est2,rep_est3, rep_est4, rep_est5, rep_est6, rep_est7,rep_est8, rep_est9)

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

bias_confusion <- lapply(rep_estimates, function(x){
  if(class(x)!= "try-error"){
    c(x[[1]] - sim[[1]]$omega)
  }
})
bias_confusion_dataframe <- do.call("rbind", bias_confusion)

#all_data <- data.frame(cbind(contrast_intercept_dataframe,bias_covariate_dataframe,bias_confusion_dataframe))
#colnames(all_data) <- c("beta21","beta31","beta32",
#                        "beta11", "beta12", "beta13",
#                        "omega11", "omega21","omega31",
#                        "omega12","omega22","omega32",
#                        "omega13","omega23", "omega33",
#                        "omega14", "omega24", "omega34")

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
 # scale_y_discrete(labels=c('beta21' = expression(Delta*(beta[0][2] - beta[0][1])),
 #                           'beta31' = expression(Delta*(beta[0][3] - beta[0][1])),
 #                           'beta32' = expression(Delta*(beta[0][3] - beta[0][2])),
 #                           'beta11'=expression(Delta*(beta[1][2] - beta[1][1])),
  #                          'beta12'=expression(Delta*(beta[1][3] - beta[1][1])),
 #                           'beta13'=expression(Delta*(beta[1][3] - beta[1][2])),
 #                           'omega11'=expression(Delta*Omega[11]),
 #                           'omega21'=expression(Delta*Omega[21]),
 #                           'omega31'=expression(Delta*Omega[31]),
 #                           'omega12'=expression(Delta*Omega[12]),
 #                           'omega22'=expression(Delta*Omega[22]),
 #                           'omega32'=expression(Delta*Omega[32]),
  #                          'omega13'=expression(Delta*Omega[13]),
 #                           'omega23'=expression(Delta*Omega[23]),
 #                           'omega33'=expression(Delta*Omega[33]),
 #                           'omega14'=expression(Delta*Omega[14]),
 #                           'omega24'=expression(Delta*Omega[24]),
 #                           'omega34'=expression(Delta*Omega[34])
  #                          ))+
  theme_bw()+
  xlab("Bias in estimates")+
  ylab("Parameters")
#print(gg1)
ggsave("bias_boxplot.png")
#z = outer(rep_estimates$beta0,rep_estimates$beta0,'-'); 
#z[lower.tri(z)];

#z1 = outer(input$fixed_effects$beta0,input$fixed_effects$beta0,'-'); 
#z1[lower.tri(z)];

#Testing the results on newdata
cov_test <- seq(-3,2.99,0.01)
cov_test <- rbind(cov_test, cov_test)
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
    beta0 = c(-3,1,2),
    beta1=c(1,0,-1),
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
library(spatstat) 

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
quantiles_beta <- apply(combine_beta, 2, function(x){quantile(x, c(0.05,0.975))})

#estimating the precision of the estimates
#combine_beta <- cbind(beta0_dataframe, beta1_dataframe)
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

#prediction of true species
n.species=3
lambda_true <- lambda_test <- prop_true <- prop_test <-prob_return_true <- prob_return_false <- prob_ret_true <- prob_ret_false <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
max_lambda_true <- max_lambda_test <- max_prob_true <- max_prob_test <- vector("numeric", length(cov_test[1,]))
omega_true = sim_test$omega
#omega_est = rep_estimates$omega
#reported_species = sim_test$Y[1,]
test_results <- function(category){
  lambda_true <- lambda_test <- prop_true <- prop_test <-prob_return_true <- prob_return_false <- prob_ret_true <- prob_ret_false <- matrix(NA, nrow=3, ncol=length(cov_test[1,]))
  max_lambda_true <- max_lambda_test <- max_prob_true <- max_prob_test <- vector("numeric", length(cov_test[1,]))
reported_species = rep(category, length(cov_test[1,]))
for(site.tag in 1:length(cov_test[1,])){
  for(spe.tag in 1:n.species){
    #True and estimated abundance
    lambda_true[spe.tag, site.tag] <- exp(input$fixed_effects$beta0[spe.tag] + input$fixed_effects$beta1[spe.tag]*cov_test[1,site.tag])
    lambda_test[spe.tag,site.tag] <- exp(median_estimates$beta0[spe.tag]+ median_estimates$beta1[spe.tag]*cov_test[1,site.tag])
    #prop_true[spe.tag, site.tag] <- proportions(lambda_true[1:n.species, site.tag])
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

g1 <- test_results(1)%>%
  dplyr:: select(cov_test, true_species1,true_species2,true_species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot()+
  geom_line(aes(x=cov_test, y=value, col=variable))+
  theme_bw()+
  xlab("")+
  ylab("Verified species prob \n from true parameters")+
  theme(legend.title = element_blank())+
  xlab("")

g2 <-ggplot(test_results(1))+
  geom_point(aes(x=cov_test, y=test_prob, col=true_prob, shape=true_prob))+
  theme_bw()+
  xlab("")+
  ylab("Predicted verified  \n species")+
  theme(legend.title = element_blank())+
  xlab("")

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

#ggplot(data, aes(x=x, y=y, group=group, fill=group)) +
#geom_density_line(stat = "identity", size=.5, alpha=0.3) +
#  scale_fill_manual(name='', values=c("1" = "green4", "2" = "red"))

gg2 <- ggarrange(g1,g3,g2, nrow=1, ncol=3, labels = c("A)","B)", "C)" ), common.legend = FALSE, legend = "none")
fig1 <- annotate_figure(gg2,
                top = text_grob("Reported = Category 1", color = "red", face = "bold", size = 14))
#ggsave("results2.png")
#ggarrange(gg1,gg2, nrow = 2, labels = c("A)"))
#ggsave("results.png")

##########
g4 <- test_results(2)%>%
  dplyr::select(cov_test, true_species1,true_species2,true_species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot()+
  geom_line(aes(x=cov_test, y=value, col=variable))+
  theme_bw()+
  xlab("")+
  ylab("Verified species prob \n from true parameters")+
  theme(legend.title = element_blank())+
  xlab("")

g5 <-ggplot(test_results(2))+
  geom_point(aes(x=cov_test, y=test_prob, col=true_prob, shape=true_prob))+
  theme_bw()+
  xlab("")+
  ylab("Predicted verified  \n species")+
  theme(legend.title = element_blank())+
  xlab("")

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

gg3 <- ggarrange(g4,g6,g5, nrow=1, ncol=3, labels = c("D)","E)", "F)" ), common.legend = FALSE, legend = "none")
fig2 <- annotate_figure(gg3,
                        top = text_grob("Reported = Category 2", color = "red", face = "bold", size = 14))
####################
g7 <- test_results(3)%>%
  dplyr:: select(cov_test, true_species1,true_species2,true_species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot()+
  geom_line(aes(x=cov_test, y=value, col=variable))+
  theme_bw()+
  xlab("")+
  ylab("Verified species prob \n from true parameters")+
  theme(legend.title = element_blank())+
  xlab("")

g8 <-ggplot(test_results(3))+
  geom_point(aes(x=cov_test, y=test_prob, col=true_prob, shape=true_prob))+
  theme_bw()+
  xlab("")+
  ylab("Predicted verified  \n species")+
  theme(legend.title = element_blank())+
  xlab("")

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

gg4 <- ggarrange(g7,g9,g8, nrow=1, ncol=3, labels = c("G)","H)", "I)" ), common.legend = FALSE, legend = "none")
fig3 <- annotate_figure(gg4,
                        top = text_grob("Reported = Category 3", color = "red", face = "bold", size = 14))

#######
g10 <- test_results(4)%>%
  dplyr::select(cov_test, true_species1,true_species2,true_species3)%>%
  reshape2::melt(id.vars=c("cov_test"), variable.names="Species")%>%
  ggplot()+
  geom_line(aes(x=cov_test, y=value, col=variable))+
  theme_bw()+
  xlab("")+
  ylab("Verified species prob \n from true parameters")+
  theme(legend.title = element_blank())+
  xlab("")

g11 <-ggplot(test_results(4))+
  geom_point(aes(x=cov_test, y=test_prob, col=true_prob, shape=true_prob))+
  theme_bw()+
  xlab("")+
  ylab("Predicted verified \n species")+
  theme(legend.title = element_blank())+
  xlab("")

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

gg5 <- ggarrange(g10,g12,g11, nrow=1, ncol=3, labels = c("J)","K)", "L)" ), common.legend = TRUE, legend = "bottom")
fig4 <- annotate_figure(gg5,
                        top = text_grob("Reported = Category 4", color = "red", face = "bold", size = 14))



library(grid)

gg_sub <- ggarrange(g3,g6,g9,g12, ncol=2, nrow=2, common.legend = TRUE, legend="bottom", labels = c("a)", "b)", "c)", "d)"))
ggpubr::annotate_figure(gg_sub, left = text_grob("Verified species probability from estimated parameters", rot = 90, vjust = 1),
                        bottom = text_grob("Covariate values"))
ggsave("results_sub.png")

gg_all <- ggarrange(fig1, fig2, fig3,fig4,
          nrow=4, ncol=1, common.legend = FALSE, legend = "bottom")
ggpubr::annotate_figure(gg_all, left = text_grob("", rot = 90, vjust = 1),
                bottom = text_grob("Covariate values"))

ggsave("results_all.png")



data1 <- data.frame(cov_test= rep(cov_test, each=3), 
                    true_est = c(lambda_true),
                    test_est = c(lambda_test),
                    species= rep(c("S1","S2","S3")), each=length(cov_test))
data1_melted <- melt(data1, id.vars = c("cov_test","species"))
vars = c(rep(c("true_S1", "true_S2", "true_S3"), each=length(cov_test)), rep(c("test_S1", "test_S2", "test_S3"), each=length(cov_test)))

data_all <- cbind(data1_melted, vars)
ggplot(data_all)+
  geom_line(aes(cov_test, value, col=as.factor(vars)))
