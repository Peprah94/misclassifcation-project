# Load the data stored from "data formatting.R" and "data_application_nimble.R"

load("data_results.RData") # results from data formatting.R
load("data_new.RData")#results from data_application_nimble.R

# packages needed to run this script
library(ggplot2)
library(pbapply)
library(purrr)
library("rnaturalearth")
library("rnaturalearthdata")
theme_set(theme_bw())
unix::rlimit_as(100*10^9)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# Number of sites 
n.sites = do.call('rbind',
                  pblapply(data_for_nimble,
                           function(x){
                             return(length(x$C))
                           }))
# Number of verified species
n.species = do.call('rbind',
                    pblapply(data_for_nimble,
                             function(x){
                               return(length(unique(c(x$C))))
                             }))

#Number of reported species
n.reported = do.call('rbind',
                     pblapply(data_for_nimble,
                              function(x){
                                return(length(unique(c(x$Y))))
                              }))

# Put all constants together
const_all <- cbind(n.sites, n.species, n.reported)

#Extract the results from the MCMC for each month
monthly_mcmc.out <- list()
for(i in 1:12){
  monthly_mcmc.out[[i]] <- rep_estimates[[i]][[9]]$summary$all.chains
}

## Classification probabilities / Confusion matrix
omega<- omega_full <- list()
for(i in 1:12){
  const <- list(n.species=const_all[i,2],n.sites=const_all[i,1], n.reported=const_all[i,3])
  omega_est <- matrix(monthly_mcmc.out[[i]][((((const$n.species*4))+ (const$n.species*(const$n.sites)))+1): ((((const$n.species*4))+ (const$n.species*(const$n.sites)))+(const$n.species*(const$n.reported))),1], nrow=const$n.species, ncol=(const$n.reported), byrow=FALSE)
  nreported =const$n.reported
  nspecies=const$n.species
  if(nreported > 5){
omega_full[[i]] = omega_est
}else{
  if(nreported == 5){
    omega_full[[i]] = cbind(omega_est,matrix(rep(0, (const$n.species)),const$n.species,1))
  }else{
    omega_full[[i]] = cbind(omega_est,matrix(rep(0, (const$n.species*(6-nreported))),const$n.species,(6-nreported)))
    }
}
  
if(nspecies > 2){
  omega[[i]] = omega_full[[i]]
}else{
  omega[[i]] = rbind(omega_full[[i]],matrix(rep(0, (6)),ncol=(6),nrow=1))
}
}
omega_res <- as.data.frame(do.call('rbind',omega))
colnames(omega_res) <- c("Danaus gilippus", "Danaus plexippus", "Lemenitis archippus", "other", "other danaus", "other lemenitis")
month <- rep(c("January", "February", "March",
      "April", "May", "June",
      "July", "August", "September",
      "October", "November", "December"), each=3)
species <- rep(c("Danaus gilippus", "Danaus plexippus", "Lemenitis archippus"),12)
misclassification_probabilities <- cbind(month, species, omega_res)
write.csv(misclassification_probabilities, "misclassification_probability.csv")

## Retrieving and saving results for the fixed effects beta0, beta1, beta2, beta3
beta0_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                if(length(x[[2]])> 2){
                                  return(x[[2]])
                                }else{
                                  return(c(x[[2]],0))
                                }
                              }))

beta1_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                if(length(x[[3]])> 2){
                                  return(x[[3]])
                                }else{
                                  return(c(x[[3]],0))
                                }

                              }))

beta2_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                if(length(x[[4]])> 2){
                                  return(x[[4]])
                                }else{
                                  return(c(x[[4]],0))
                                }
                              }))

beta3_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                if(length(x[[5]])> 2){
                                  return(x[[5]])
                                }else{
                                  return(c(x[[5]],0))
                                }
                              }))
month <- c("January", "February", "March",
               "April", "May", "June",
               "July", "August", "September",
               "October", "November", "December")
betas_res <- cbind(month,beta0_res,
                   beta1_res,
                   beta2_res,
                   beta3_res)
write.csv(betas_res, "betas.csv")

# Function that estimates the proportion of the entire dataset.
proportion_estimate <- function(data){
  ret <- matrix(NA, nrow=nrow(data), ncol = ncol(data))
  for(i in 1:nrow(data)){
    ret[i,] <- proportions(data[i,])
  }
  return(ret)
}

###########################################
# Prediction for 2019 data that needed ID known as validation data
##########################################
results <- max_species <- list()
for(month in 1:12){
  #for(i in 1:n.species){
  print(month)
  print(rowSums(omega[[month]]))
  if(length(data_for_nimble[[month]]$validation_C) >0){
  pred1 <- exp(beta0_res[month,1]+ beta1_res[month,1]%*%data_for_nimble[[month]]$temp_validation + beta2_res[month,1]%*%data_for_nimble[[month]]$altitude_validation+ beta3_res[month,1]%*%data_for_nimble[[month]]$validation_agree)
  pred2 <- exp(beta0_res[month,2]+ beta1_res[month,2]%*%data_for_nimble[[month]]$temp_validation + beta2_res[month,2]%*%data_for_nimble[[month]]$altitude_validation+ beta3_res[month,2]%*%data_for_nimble[[month]]$validation_agree)
  pred3 <- exp(beta0_res[month,3]+ beta1_res[month,3]%*%data_for_nimble[[month]]$temp_validation + beta2_res[month,3]%*%data_for_nimble[[month]]$altitude_validation+ beta3_res[month,3]%*%data_for_nimble[[month]]$validation_agree)
  }else{
    pred1 <- pred2 <- pred3 <- 1
    
  }
  pred_matrix <- proportion_estimate(cbind(t(pred1), t(pred2), t(pred3)))
  pred_prop_matrix <- pred_matrix/rowSums(pred_matrix)
  
  res <- list()
  if(length(data_for_nimble[[month]]$validation_C) >0){
  for(i in 1:(length(data_for_nimble[[month]]$validation_Y))){
    print(i)
      coords_validation = data_for_nimble[[month]]$coords_validation
      res1 <-  pred_prop_matrix[i,]* omega[[month]][,(data_for_nimble[[month]]$validation_Y[i])]/sum(pred_prop_matrix[i,]* omega[[month]][,(data_for_nimble[[month]]$validation_Y[i])])
      max_species <- which.max(res1)
      res[[i]] <- c(month,
                    c(res1),
                    data_for_nimble[[month]]$validation_Y[i], 
                    data_for_nimble[[month]]$validation_C[i],
                    max_species,
                    coords_validation[i,]$X1,
                    coords_validation[i,]$X2)
  }
      }else{
        for(i in 1:2){ #create a dummy iteration for flatten to work
     coords_validation =data.frame(X1 =0, X2=0)
     res1 <- rep(0,3)
     res[[i]] <- c(month,
                   c(res1),
                   0, 
                   0,
                  99,
                   coords_validation$X1,
                   coords_validation$X2)
     }
    }
  results[[month]] <- res
}
pred_probs = flatten(results)
full_data <- data.frame(do.call("rbind", pred_probs))
colnames(full_data) <- c("month", "Danaus gilipus", "Danaus plexippus", "Limenitis archippus", "Reported Species","truth", "max species","longitude", "latitude")
sum(full_data$truth==full_data$`max species`)/nrow(full_data) #proportion of max species that correspond to the truth

#save results
write.csv(full_data, "results_from_data.csv")

#monthly percentage of correct identifications
validation_prob <- full_data %>%
  select(truth, `max species`, month) %>%
  dplyr:: group_by(month) %>%
  summarise(probability = mean(truth==`max species`, na.rm=TRUE) )

#Plot of results
melted_data <- reshape2::melt(full_data, id.vars=c("month","Reported Species","truth","max species","longitude", "latitude"))
ggplot(data = world)+
  geom_sf()+
  geom_point(data =melted_data, aes(as.numeric(longitude), as.numeric(latitude), col=value))+
  coord_sf(xlim = c(-125, -68), ylim = c(17.5, 46.9), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")+
  facet_wrap(vars(variable))

ggplot(data = world)+
  geom_sf()+
  geom_point(data =melted_data, aes(as.numeric(longitude), as.numeric(latitude), col=as.factor(truth)))+
  coord_sf(xlim = c(-125, -68), ylim = c(17.5, 46.9), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")#+
  #facet_wrap(~month, ncol=3)+
  #theme(legend.position = "bottom")

###########################################
# Prediction for 2020 research grade data
##########################################
results_20 <- list()
for(month in 1:12){
  #for(i in 1:n.species){
  print(month)
  pred1 <- exp(beta0_res[month,1]+ beta1_res[month,1]%*%data_for_nimble[[month]]$temp_20 + beta2_res[month,1]%*%data_for_nimble[[month]]$altitude_20+ beta3_res[month,1]%*%data_for_nimble[[month]]$research20_agree)
  pred2 <- exp(beta0_res[month,2]+ beta1_res[month,2]%*%data_for_nimble[[month]]$temp_20 + beta2_res[month,2]%*%data_for_nimble[[month]]$altitude_20+ beta3_res[month,2]%*%data_for_nimble[[month]]$research20_agree)
  pred3 <- exp(beta0_res[month,3]+ beta1_res[month,3]%*%data_for_nimble[[month]]$temp_20 + beta2_res[month,3]%*%data_for_nimble[[month]]$altitude_20+ beta3_res[month,3]%*%data_for_nimble[[month]]$research20_agree)
  dat <- cbind(t(pred1), t(pred2), t(pred3))
  pred_matrix <- proportion_estimate(dat)
  pred_prop_matrix <- pred_matrix/rowSums(pred_matrix)
  res <- list()
  if(length(data_for_nimble[[month]]$research_C20) >0){
    for(i in 1:(length(data_for_nimble[[month]]$research_Y20))){
      print(i)
      #iter <- month*i
      #print(iter)
      coords_validation = data_for_nimble[[month]]$coords_20
      res1 <-  pred_prop_matrix[i,]* omega[[month]][,(data_for_nimble[[month]]$research_Y20[i])]/sum(pred_prop_matrix[i,]* omega[[month]][,(data_for_nimble[[month]]$research_Y20[i])])
      max_species <- which.max(res1)
      res[[i]] <- c(month,
                    c(res1),
                    data_for_nimble[[month]]$research_Y20[i], 
                    data_for_nimble[[month]]$research_C20[i],
                    max_species,
                    coords_validation[i,]$X1,
                    coords_validation[i,]$X2)
    }
  }else{
    for(i in 1:2){ #create a dummy iteration for flatten to work
      coords_validation =data.frame(X1 =0, X2=0)
      res1 <- rep(0,3)
      res[[i]] <- c(month,
                    c(res1),
                    0, 
                    0,
                    99,
                    coords_validation$X1,
                    coords_validation$X2)
    }
  }
  results_20[[month]] <- res
}

pred_probs = flatten(results_20)
full_data <- data.frame(do.call("rbind", pred_probs))
colnames(full_data) <- c("month", "Danaus gilipus", "Danaus plexippus", "Limenitis archippus", "truth","Reported Species", "max_species","longitude", "latitude")
#percentage of correct identifications
sum(full_data$truth==full_data$max_species)/nrow(full_data)

#monthly percentage of correct identifications
research_2020 <-full_data %>%
  select(truth, max_species, month) %>%
 dplyr:: group_by(month) %>%
  summarise(probability = mean(truth==max_species, na.rm=TRUE) )

#save the data
write.csv(full_data, "results_from_data_2020.csv")

#Plot the results
melted_data <- reshape2::melt(full_data, id.vars=c("month","Reported Species","truth","max_species","longitude", "latitude"))
gg1 <- list()
for(i in 1:12){
  melted_data_new <- melted_data %>%
    filter(month==i)
 gg1[[i]] <-  ggplot(data = world)+
  geom_sf()+
  geom_point(data =melted_data_new, aes(as.numeric(longitude), as.numeric(latitude), color=value))+
  coord_sf(xlim = c(-125, -68), ylim = c(17.5, 46.9), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")+
   labs(title=paste("Month ",i) ,color="Predicted probability\n")+
  facet_wrap(~variable, ncol=3)+
  theme(legend.position = "bottom")
}
ggarrange(gg1[[1]],gg1[[2]], gg1[[3]], gg1[[4]],
          gg1[[5]],gg1[[6]],gg1[[7]],gg1[[8]],
          gg1[[9]],gg1[[10]],gg1[[11]],gg1[[12]],
          common.legend = TRUE, nrow=3, legend="bottom")
 ggplot(data = world)+
  geom_sf()+
  geom_point(data =melted_data, aes(as.numeric(longitude), as.numeric(latitude), col=as.factor(`Reported Species`)))+
  coord_sf(xlim = c(-125, -68), ylim = c(17.5, 46.9), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")+
  facet_wrap(~month, ncol  =6)+
  theme(legend.position = "bottom")
 
 
 #Combine all probabilities and plot them
 all_probability <- rbind(validation_prob ,research_2020)
groups <- as.factor(rep(c("needs ID 2019", "research grade 2020"), each=12))
all_data <- cbind(groups,all_probability)
ggplot(data = all_data, aes(x=as.factor(month), y=probability,color=as.factor(groups), group=as.factor(groups)))+
  geom_point()+
  geom_line()+
  xlab("Month")+
  ylab("Probability of correct classification")+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.2))


            
