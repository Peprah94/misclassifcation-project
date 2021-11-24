load("/Volumes/kwakupa/misclassification/data_results.RData")
load("~/Documents/GitHub/Thesis/Misclassification/Misclassification/data_new.RData")
library(ggplot2)
library(pbapply)
library(purrr)
library("rnaturalearth")
library("rnaturalearthdata")
theme_set(theme_bw())
unix::rlimit_as(100*10^9)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

const <- list(n.sites=1303,
              n.species=3,
              #n.sites = dim_data[3], 
              #n.species = (dim_data[1]), 
              n.reported=4)

beta0_res <- do.call('rbind',
        pblapply(rep_estimates, 
                 function(x){
                   return(x[[2]])
                 }))

beta1_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                return(x[[3]])
                              }))

beta2_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                return(x[[4]])
                              }))

beta3_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                return(x[[5]])
                              }))

betas_res <- cbind(beta0_res,
                   beta1_res,
                   beta2_res,
                   beta3_res)

omega_res <- do.call('rbind',
                     pblapply(rep_estimates, 
                              function(x){
                                return(x[[1]])}))
 
#Does not sum to 1, check out why that is happening in the data                                                             }))
rowSums(omega_res)

# Prediction of lambda
#pred_matrix <- pred_prop_matrix <- gamma <- list()
results <- list()
for(month in 1:12){
  #for(i in 1:n.species){
  print(month)
    pred1 <- exp(beta0_res[month,1]+ beta1_res[month,1]%*%data_for_nimble[[month]]$temp_validation + beta2_res[month,1]%*%data_for_nimble[[month]]$altitude_validation+ beta3_res[month,1]%*%data_for_nimble[[month]]$validation_agree)
    pred2 <- exp(beta0_res[month,2]+ beta1_res[month,2]%*%data_for_nimble[[month]]$temp_validation + beta2_res[month,2]%*%data_for_nimble[[month]]$altitude_validation+ beta3_res[month,2]%*%data_for_nimble[[month]]$validation_agree)
    pred3 <- exp(beta0_res[month,3]+ beta1_res[month,3]%*%data_for_nimble[[month]]$temp_validation + beta2_res[month,3]%*%data_for_nimble[[month]]$altitude_validation+ beta3_res[month,3]%*%data_for_nimble[[month]]$validation_agree)
    pred_matrix <- proportions(cbind(t(pred1), t(pred2), t(pred3)))
    pred_prop_matrix <- pred_matrix/rowSums(pred_matrix)
    res <- list()
    for(i in 1:(length(data_for_nimble[[month]]$validation_Y))){
   print(i)
      #iter <- month*i
      #print(iter)
      res1 <-  pred_prop_matrix[i,]* rep_estimates[[month]][[1]][,(data_for_nimble[[month]]$validation_Y[i])]/sum(pred_prop_matrix[i,]* rep_estimates[[month]][[1]][,(data_for_nimble[[month]]$validation_Y[i])])
    res[[i]] <- c(month,
                  c(res1),
                  data_for_nimble[[month]]$validation_Y[i], 
                  data_for_nimble[[month]]$coords_validation[i,]$X1,
                  data_for_nimble[[month]]$coords_validation[i,]$X2
                                                                              )
      }
    results[[month]] <- res
}

results_19 <- list()
for(month in 1:12){
  #for(i in 1:n.species){
  print(month)
  pred1 <- exp(beta0_res[month,1]+ beta1_res[month,1]%*%data_for_nimble[[month]]$temp_19 + beta2_res[month,1]%*%data_for_nimble[[month]]$altitude_19+ beta3_res[month,1]%*%data_for_nimble[[month]]$research19_agree)
  pred2 <- exp(beta0_res[month,2]+ beta1_res[month,2]%*%data_for_nimble[[month]]$temp_19 + beta2_res[month,2]%*%data_for_nimble[[month]]$altitude_19+ beta3_res[month,2]%*%data_for_nimble[[month]]$research19_agree)
  pred3 <- exp(beta0_res[month,3]+ beta1_res[month,3]%*%data_for_nimble[[month]]$temp_19 + beta2_res[month,3]%*%data_for_nimble[[month]]$altitude_19+ beta3_res[month,3]%*%data_for_nimble[[month]]$research19_agree)
  pred_matrix <- proportions(cbind(t(pred1), t(pred2), t(pred3)))
  pred_prop_matrix <- pred_matrix/rowSums(pred_matrix)
  res <- list()
  for(i in 1:(length(data_for_nimble[[month]]$research_Y19))){
    print(i)
    #iter <- month*i
    #print(iter)
    res1 <-  pred_prop_matrix[i,]* rep_estimates[[month]][[1]][,(data_for_nimble[[month]]$research_Y19[i])]/sum(pred_prop_matrix[i,]* rep_estimates[[month]][[1]][,(data_for_nimble[[month]]$research_Y19[i])])
    res[[i]] <- c(month,
                  c(res1),
                  data_for_nimble[[month]]$research_Y19[i], 
                  data_for_nimble[[month]]$coords_19[i,]$X1,
                  data_for_nimble[[month]]$coords_19[i,]$X2
    )
  }
  results_19[[month]] <- res
}


pred_probs = flatten(results)
full_data <- data.frame(do.call("rbind", pred_probs))
colnames(full_data) <- c("month", "species 1", "Species 2", "Species 3", "Reported Species", "longitude", "latitude")
write.csv(full_data, "results_from_data.csv")

melted_data <- reshape2::melt(full_data, id.vars=c("month","Reported Species","longitude", "latitude"))
ggplot(data = world)+
  geom_sf()+
  geom_point(data =melted_data, aes(as.numeric(longitude), as.numeric(latitude), col=value))+
  coord_sf(xlim = c(-125, -68), ylim = c(17.5, 46.9), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")+
  facet_wrap(vars(variable))
