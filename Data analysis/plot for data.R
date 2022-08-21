setwd("/Volumes/kwakupa/misclassification/new")
load("/Volumes/kwakupa/misclassification/new/nimble_data.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_constant1.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_intercept1.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_main1.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_only_principal_cov1.RData")
only_principal_cov <- variable
load("/Volumes/kwakupa/misclassification/new/estimated_data_variable1.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_op_cov_finter1.RData")

all_data <- list(variable, 
                 constant, 
                 intercept,
                 main, 
                 op_cov_finter,
                 only_principal_cov)
#packages
library(ggplot2)
library(dplyr)
theme_set(theme_bw())

#plot_fnx <- function(gamma0, user_effect, N){
#  alpha <- exp(gamma0 + user_effect[N])
#  omega <- proportions(alpha, margin = 1)
#  ret <- c(omega[1,1], omega[2,2], omega[3,3])
#  return(ret)
#}


plot_fnx <- function(gamma0, user_effect, N){
  alpha <- omega <-  array(NA, dim = c(3,6, length(N)))
  ret <- matrix(NA, nrow = length(N), ncol = 3)
 # if(type == "annex"){

  for(site.tag in 1:length(N)){
    for(spe.tag in 1:3){
  alpha[spe.tag,6 , site.tag] <- 1
  alpha[spe.tag,1:5 , site.tag] <- exp(gamma0[spe.tag, ] + user_effect[N])
  omega[spe.tag, , site.tag] <- proportions(alpha[spe.tag, , site.tag])
    }
    ret[site.tag, ] <- c(omega[1,1, site.tag], omega[2,2, site.tag], omega[3,3, site.tag])
  }
  #}else{
    #for(site.tag in 1:length(N)){
     # for(spe.tag in 1:3){
      #  alpha[spe.tag,6 , site.tag] <- 1
     #   alpha[spe.tag,1:5 , site.tag] <- exp(gamma0[spe.tag, ] + gamma1[spe.tag, ] * log(N))
      #  omega[spe.tag, , site.tag] <- proportions(alpha[spe.tag, , site.tag])
     # }
     # ret[site.tag, ] <- c(omega[1,1, site.tag], omega[2,2, site.tag], omega[3,3, site.tag])
    #}
  #}
  return(ret)
}

subset_parameters_annex <- function(data, subset_index){
  ret <-  lapply(data, function(x){
    d <- (x[[subset_index]])
  })
  return(ret)
}

subset_parameters <- function(data, subset_index){
  ret <-  lapply(data, function(x){
    d <- c(x[[subset_index]])
  })
  return(ret)
}

gamma0 <- subset_parameters_annex(all_data, 7)
gamma1 <- subset_parameters_annex(all_data, 8)
user_effect <- subset_parameters_annex(all_data, 12)

omega <- lapply(all_data, function(x){
  #i <- 1
  subset_index <- c(1,5,9, 19, 23, 27,
                    37,41, 45, 55,59, 63,
                    73,77,81, 91, 95, 99)
  
  tt <- sum(!is.na(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 2][subset_index]))
  if(tt ==3){
  main<- rep(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 2][c(1,5,9)],6)
 lower <- rep(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 4][c(1,5,9)],6)
 upper <- rep(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 5][c(1,5,9)],6)
  }else{
    main<- x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 2][subset_index]
    lower <- x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 4][subset_index]
    upper <- x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 5][subset_index] 
  }
 N <- rep(1:6, each = 3)
 #i  <- i +1
 #print(i)
 ret <- cbind(N, main, lower, upper)
 return(ret)
})%>%
  do.call("rbind", .)%>%
  data.frame()%>%
  dplyr::mutate(Method = c(rep("CCM", 18),
                    rep("COM", 18),
                    rep("CIM", 18),
                    rep("MCM", 18),
                    rep("CSICM", 18),
                    rep("CSCM", 18)),
                Species = rep(c("Queen","Monarch", "viceroy"),36))%>%
  ggplot()+
  #geom_point(aes(x = N, y = value, col = Method))+
  geom_line(aes(x = as.factor(N), y = main, col = Method, group = Method))+
  geom_ribbon(aes(x = as.factor(N), ymin= lower, ymax= upper, fill= Method, group = Method), alpha = 0.2)+
  ylab("Probability of correct \n classification")+
  xlab("experience")+
  #scale_x_discrete("N", labels = as.character(N), breaks = N)
  facet_wrap(~Species)
omega
experience <- ggplot(data = data.frame(experience = dataout$experience))+
  geom_bar(aes(x = as.factor(experience)))+
  xlab("experience")

fig <- ggpubr::ggarrange(omega, experience, nrow =2 , ncol = 1)

ggsave("omega.png", plot = fig, dpi = 320, 
       height = 10, width = 18, units = "cm")







#all_data <- list(variable, intercept, constant, main)

accuracy <- subset_parameters(all_data, 9)%>%
  unlist()

precision <- subset_parameters(all_data, 10)%>%
  unlist()

recall <- subset_parameters(all_data, 13)%>%
  unlist()

method = c("CCM", "COM", "CIM", "MCM","CSICM", "CSCM")

pred_performance <- cbind(method, accuracy, precision, recall)
write.csv(pred_performance,"pred_performance.csv", row.names = FALSE)

p <- subset_parameters(all_data, 11)%>%
  unlist()

beta0 <- subset_parameters(all_data, 1)%>%
  unlist()

beta1 <- subset_parameters(all_data, 2)%>%
  unlist()

beta2 <- subset_parameters(all_data, 3)%>%
  unlist()
beta3 <- subset_parameters(all_data, 4)%>%
  unlist()
beta4 <- subset_parameters(all_data, 5)%>%
  unlist()
beta5 <- subset_parameters(all_data, 6)%>%
  unlist()

parameters_estimates <- cbind(beta0, beta1, beta2, beta3, beta4, beta5)%>%
  data.frame()%>%
 dplyr::mutate(method = c(rep("CCM", 2),
                    rep("COM", 2),
                    rep("CIM", 2),
                    rep("MCM", 2),
                    rep("CSICM", 2),
                    rep("CSCM",2)),
               species = rep(c("Monarch", "viceroy"), 6))%>%
  reshape2::melt(id.vars = c("method", "species"))%>%
  reshape2::dcast(., species+variable ~ method)%>%
  dplyr::arrange(variable)%>%
  dplyr::mutate(parameter = rep(c("Intercept", "Altitude",
                                  "Distance to road", "Precipitation",
                                  "Mean temperature", "Experience"), each = 2))
write.csv(parameters_estimates,"parameters_estimates.csv", row.names = FALSE)



###############
#   Classification process paramaters
##############
gamma0 <- subset_parameters_annex(all_data, 7)%>%
  unlist()
#gamma1 <- subset_parameters_annex(all_data, 8)
user_effect <- subset_parameters_annex(all_data, 12)%>%
  unlist()

misclass_pars <- data.frame(value = gamma0)%>%
  dplyr::mutate(method = c(rep("CCM", 15),
                           rep("COM", 15),
                           rep("CIM", 15),
                           rep("MCM", 15),
                           rep("CSICM", 15),
                           rep("CSCM", 15)),
                verified_species = rep(c("Queen","Monarch", "viceroy"),30),
                classified_states = rep(rep(c("Queen","Monarch", "viceroy",
                                              "other", "other_danaus"), each = 3), 6))%>%
  reshape2::dcast(., verified_species + classified_states ~ method)
write.csv(misclass_pars,"misclassification_intercept.csv", row.names = FALSE)





user_effect <- data.frame(value = user_effect)%>%
  dplyr::mutate(method = c(rep("CCM",6),
                           rep("COM", 6),
                           rep("CIM", 6),
                           rep("MCM", 6),
                           rep("CSIFM", 6)))
