#load("/Volumes/kwakupa/testrun/estimateddata_new.RData")
#sim <- load("/Volumes/kwakupa/testrun/simulated data.RData")
n.sites=900; n.visit=1; n.species=3
results <- rep_estimates
x0= sim[[1]]$locs[,1]
y0 = sim[[1]]$locs[,2]
gridcov = sim[[1]]$cov
logLam1 = as.numeric(log(results$lambda[1,]))
logLam2 = as.numeric(log(results$lambda[2,]))
logLam3 = as.numeric(log(results$lambda[3,]))
library(ggplot2)
library(ggpubr)

data1 <- data.frame(x=x0, y=y0, loglambda1=logLam1, loglambda2 = logLam2, loglambda3=logLam3)
g1 <- ggplot(data1, aes(x0,y0))+
  geom_raster(aes(fill=loglambda1))+
  scale_fill_continuous(breaks=c(0,2,4,6))+
  scale_fill_gradientn(colors = terrain.colors(20))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("A) Predicted Log Intensity of S1")

g2 <- ggplot(data1, aes(x0,y0))+
  geom_raster(aes(fill=loglambda2))+
  scale_fill_continuous(breaks=c(0,2,4,6))+
  scale_fill_gradientn(colors = terrain.colors(20))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("B)Predicted Log Intensity of S2")

g3 <- ggplot(data1, aes(x0,y0))+
  geom_raster(aes(fill=loglambda3))+
  scale_fill_continuous(breaks=c(0,2,4,6))+
  scale_fill_gradientn(colors = terrain.colors(20))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("C)Predicted Log Intensity of S3")
ggarrange(g1,g2,g3, nrow=2, ncol=2)

###################
# predicting true species
#####################

predicted_omega <- rep_estimates$omega
predicted_prop <- rep_estimates$prop

prob_for_true <- function(predicted_omega, predicted_prop, sim,val){
true_pred_prob <- matrix(NA, nrow=3, ncol=n.sites)
index_val <- c()
for(i in 1:n.species){
for(j in 1:n.sites){
  if(!is.na(sim[[1]]$Y[j])){
  index_val[j] <- sim[[1]]$Y[j]}
  else{
    index_val[j] <- val
  }
  true_pred_prob[i,j] <- (predicted_prop[i,j]*predicted_omega[i,index_val[j]])
}
}
true_pred_prob1 <- t(t(true_pred_prob)/rowSums(t(true_pred_prob)))
#colSums(t(true_pred_prob1))
return(true_pred_prob1)
}

max_prob_true <- function(predicted_omega, predicted_prop, sim,val){
  max_prob <- c()
  prob <- prob_for_true(predicted_omega, predicted_prop, sim,val)
for(i in 1:n.sites){
  max_prob[i]<- which.max(prob[,i])
}
  return(max_prob)
}

#prob_for_true(predicted_omega, predicted_prop, sim, 1)
index1 <- as.factor(max_prob_true(predicted_omega, predicted_prop, sim, 1))
index2 <-as.factor(max_prob_true(predicted_omega, predicted_prop, sim, 2))
index3 <- as.factor(max_prob_true(predicted_omega, predicted_prop, sim, 3))
index4 <- as.factor(max_prob_true(predicted_omega, predicted_prop, sim, 4))

data2 <- data.frame(x=x0, y=y0, index1=index1, index2 = index2, index3=index3, index4=index4)
g1 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index1))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("A) Predicted VS using \n maxprob (NAs=1)")

g2 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index2))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("C) Predicted VS using \n maxprob (NAs=2)")

g3 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index3))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("E) Predicted VS using \n maxprob (NAs=3)")

g4 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index4))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("G) Predicted VS using \n maxprob (NAs=4)")

#ggarrange(g1,g2,g3,g4, nrow=2, ncol=2)


#######
sim_true <- function(predicted_omega, predicted_prop, sim,val){
C <- array(NA, dim=c(n.species, n.visit,n.sites)) # Obs detect/non-detect
numb <- array(NA, dim=c(n.visit, n.sites))
prop <- prob_for_true(predicted_omega, predicted_prop, sim,val)
for(site.tag in 1:n.sites){
  for(visit.tag in 1:n.visit){
    C[,visit.tag, site.tag] <- rmultinom(1,1, prop[1:n.species,site.tag])
  }
}


# Which species were present
for(visit.tag in 1:n.visit){
  for(site.tag in 1:n.sites){
    numb[visit.tag,site.tag] <-which(C[1:n.species,visit.tag,site.tag]==1)
  }
}
return(numb)
}
index1 <- as.factor(sim_true(predicted_omega, predicted_prop, sim, 1))
index2 <-as.factor(sim_true(predicted_omega, predicted_prop, sim, 2))
index3 <- as.factor(sim_true(predicted_omega, predicted_prop, sim, 3))
index4 <- as.factor(sim_true(predicted_omega, predicted_prop, sim, 4))


data2 <- data.frame(x=x0, y=y0, index1=index1, index2 = index2, index3=index3, index4=index4)
g5 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index1))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("B) Predicted VS using \n simulation (NAs=1)")

g6 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index2))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("D) Predicted VS using \n simulation (NAs=2)")

g7 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index3))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("F) Predicted VS using \n simulation (NAs=3)")

g8 <- ggplot(data2, aes(x0,y0))+
  geom_raster(aes(fill=index4))+
  scale_fill_manual(values = c("green", "red", "blue"),labels=c("S1", "S2", "S3") ,guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("H) Predicted VS using \n simulation (NAs=4)")

ggarrange(g1,g5,g2,g6,g3,g7,g4,g8, nrow=2, ncol=4, common.legend = TRUE, legend = "bottom")

