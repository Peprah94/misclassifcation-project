# Function to generate data for simulations

# Function returns a simulated log cox gaussian process with or without the mean of this process related to an environmental variable. 

## Inputs:

#dim - dimensions of the window
#lambda - mean of random field
#env.beta - coefficient of env variable
#plotdat - logical, plot data or not
#seed - value to set seed
#sigma2x and kappa = parameters for the matern field

## Extra functionality

# - env correlated with spatial bias or not


### SImulation of covariates
set.seed(1994)
library(spatstat) 
# owin creates an object of class "owin" which is an observation window in 2D
# specify x and y coordinates
#win <-owin(c(-1,4), c(-1,4))
dim = c(20,20)
win <- owin(c(0,dim[1]), c(0,dim[2])) # keeping this rectangular is really important to check for errors in the code - otherwise easy to get x and y confused
# set number of pixels
spatstat.options(npixel=c(dim[1],dim[2]))

## Creating environmental covariate and altering truth based on this
# using Simpson tutorial again

# can add covariate at the beginning - create arificial covariate
# create y0 and x0 separately as rectangle
y0 <- seq(win$yrange[1], win$yrange[2],
          length=spatstat.options()$npixel[2])
x0 <- seq(win$xrange[1], win$xrange[2],
          length=spatstat.options()$npixel[1])
# bit of a fudge but rounding this does work to give 3 levels (0,1,2)
#gridcov <- round(outer(y0, x0, function(x,y) cos(x/6) - sin((y/6)-2))) 
#x0 <- seq(-1, 4, length = 100)
#y0 <- seq(-1,4, length = 100)
locs <- expand.grid(x0,y0)
# to make this more complex, changed to a continuous gradient from the bottom to the top with values from 0 to 1
multiplier <- 1/dim[2]

gridcov1 <- outer(y0,x0, function (x,y) multiplier*x + 0*y)
#gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
#gridcov <- rnorm(n.sites, 0,1)
#gridcov1 <- rnorm(n.sites,0,1)
gridcov <- runif(dim[1]*dim[2], -2,2)
cov <- rbind(c(t(gridcov)), c(t(gridcov1)))

###Set parameters

#genData
library(pbapply)

input <- list(
  constants=list(
    n.species=3,
    #n.sites=nrow(locs),
    n.visit=1, 
    p.tag= c(0.7,0.8,0.90),
    dim = c(20,20),
  sample.size = 400
),
plotdata = list(
  log_lambda=FALSE,
  prop_plot=FALSE),
fixed_effects = list(
  beta0 = c(2,3,2),
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


genData <- function(input, seed=1){ # made all of the inputs null so can see everything feeds through
  dim= input$constants$dim
  n.sites = dim[1]*dim[2]
  n.visit=input$constants$n.visit
  p.tag = input$constants$p.tag
  n.species=input$constants$n.species
  library(spatstat) 
  # owin creates an object of class "owin" which is an observation window in 2D
  # specify x and y coordinates
  #win <-owin(c(-1,4), c(-1,4))
  win <- owin(c(0,dim[1]), c(0,dim[2])) # keeping this rectangular is really important to check for errors in the code - otherwise easy to get x and y confused
  # set number of pixels
  spatstat.options(npixel=c(dim[1],dim[2]))
  
  ## Creating environmental covariate and altering truth based on this
  # using Simpson tutorial again
  
  # can add covariate at the beginning - create arificial covariate
  # create y0 and x0 separately as rectangle
 # y0 <- seq(win$yrange[1], win$yrange[2],
          #  length=spatstat.options()$npixel[2])
  #x0 <- seq(win$xrange[1], win$xrange[2],
          #  length=spatstat.options()$npixel[1])
  # bit of a fudge but rounding this does work to give 3 levels (0,1,2)
  #gridcov <- round(outer(y0, x0, function(x,y) cos(x/5) - sin((y/5)-2))) 
  #x0 <- seq(-1, 4, length = 100)
  #y0 <- seq(-1,4, length = 100)
  locs <- expand.grid(x0,y0)
  # to make this more complex, changed to a continuous gradient from the bottom to the top with values from 0 to 1
  #multiplier <- 1/dim[2]
  
  #gridcov <- outer(y0,x0, function (x,y) multiplier*x + 0*y)
  #gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
  #gridcov <- rnorm(n.sites, 0,1)
 # gridcov1 <- rnorm(n.sites,0,1)
  
  
  beta0 <- input$fixed_effects$beta0 # intercept/mu (increased to 5 to increase no.obs)
  beta1 <- input$fixed_effects$beta1 # slope of relationship to environment
  beta2 <- input$fixed_effects$beta2
  

  cov <- input$covariates$cov
  lambda <- matrix(NA, nrow=n.species, ncol=n.sites)
  for(spe.tag in 1:n.species){
    for(site.tag in 1:n.sites){
      lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag]*cov[1,site.tag])#+beta2[spe.tag]*cov[2,site.tag])
    }
  }
  
  
  if (input$plotdata$log_lambda==TRUE){
    library(ggplot2)
    library(ggpubr)
    
    #data <- data.frame(x=locs$Var1, y=locs$Var2, value= c(t(gridcov)), value1 = c(t(gridcov1)))
   # g1 <- ggplot(data, aes(x,y))+
   #   geom_raster(aes(fill=value))+
    #  scale_fill_gradientn(colors = terrain.colors(20))+
    #  theme_bw()+
    #  coord_fixed()+
    #  theme(legend.title = element_blank())+
    #  xlab("")+
    #  ylab("")+
    #  ggtitle("A) Covariate")
    
    #g5 <- ggplot(data, aes(x,y))+
     # geom_raster(aes(fill=value1))+
      #scale_fill_gradientn(colors = terrain.colors(20))+
      #theme_bw()+
     # coord_fixed()+
     # theme(legend.title = element_blank())+
     # xlab("")+
    #  ylab("")+
    #  ggtitle("A) Covariate")
    
    data1 <- data.frame(x=locs$Var1, y=locs$Var2, loglambda1=lambda[1,], loglambda2 = lambda[2,], loglambda3=lambda[3,])
    g2 <- ggplot(data1, aes(x,y))+
      geom_raster(aes(fill=loglambda1))+
      scale_fill_continuous(breaks=c(0,2,4,6))+
      scale_fill_gradientn(colors = terrain.colors(20))+
      theme_bw()+
      coord_fixed()+
      theme(legend.title = element_blank())+
      xlab("")+
      ylab("")+
      ggtitle("B) Log Intensity of S1")
    
    g3 <- ggplot(data1, aes(x,y))+
      geom_raster(aes(fill=loglambda2))+
      scale_fill_continuous(breaks=c(0,2,4,6))+
      scale_fill_gradientn(colors = terrain.colors(20))+
      theme_bw()+
      coord_fixed()+
      theme(legend.title = element_blank())+
      xlab("")+
      ylab("")+
      ggtitle("C) Log Intensity of S2")
    
    g4 <- ggplot(data1, aes(x,y))+
      geom_raster(aes(fill=loglambda3))+
      scale_fill_continuous(breaks=c(0,2,4,6))+
      scale_fill_gradientn(colors = terrain.colors(20))+
      theme_bw()+
      coord_fixed()+
      theme(legend.title = element_blank())+
      xlab("")+
      ylab("")+
      ggtitle("D) Log Intensity of S3")
    
   print(ggarrange(g1,g2,g3,g4,g5, nrow=2, ncol=3, common.legend = FALSE, legend="bottom"))
    }
   # points(xy, pch=19)

  
  
 #Combine lambda for all the species
com_lambda <- prop <-  matrix(NA, ncol=n.sites ,nrow=n.species) 
#for(spe.tag in 1:n.species){
#  com_lambda[spe.tag,] <- c(t(rf.s[[spe.tag]]))
#}
  max_prob = c()
for(site.tag in 1:n.sites){
  prop[,site.tag]<- proportions(lambda[1:n.species, site.tag])
  max_prob[site.tag] <- which.max(prop[, site.tag])
}

#max_prob <- matrix(max_prob, nrow= dim[1], ncol=dim[2], byrow=F)


  
#Check if the probabilities are correct and sum up to 1
for(site.tag in n.sites){
  for(spe.tag in 1:n.species){
    if(any(prop[spe.tag, site.tag] <0) | any(prop[spe.tag, site.tag] >1) ) stop("Entry in Probability wrong")
    if(round(sum(prop[1:n.species,site.tag ]),0) != 1) stop("Sum of probabilities is not equal to 1")
  }
}

z <- array(NA, dim=c(n.species, n.visit, n.sites)) #Detect/Non-Detection
C <- array(NA, dim=c(n.visit,n.sites)) # Obs detect/non-detect
Y <- numb <- numb_with_NA <- array(NA, dim=c(n.visit, n.sites)) # Verified species with Machine Learning
omega <- array(NA, dim= c(n.species, n.species+1)) # Classification Confusion matrix

for(site.tag in 1:n.sites){
  #   for(visit.tag in 1:n.visit){
  omega[1,2] <- 0.05
  omega[1,3] <- 0.13
  omega[2,1] <- 0.10
  omega[2,3] <- 0.05
  omega[3,1] <- 0
  omega[3,2] <- 0
  for(spe.tag in 1:n.species){
    omega[spe.tag, spe.tag] <- p.tag[spe.tag]
    omega[spe.tag, n.species+1] <- 1 - sum(omega[spe.tag, 1:n.species])  
    #}
  }
}

# Check if the rowSums are 1
for(spe.tag in 1:n.species){
  for(visit.tag in 1:n.visit){
    for(site.tag in n.sites){
      if(sum(omega[spe.tag,]) != 1) stop("Sum of confusion matrix is not equal to 1")
      if(any(omega[spe.tag,] <0) | any(omega[spe.tag,] >1) ) stop("Entry with Confusion matrix")
    }
  }
}

#Verified species
library(boot)
library(extraDistr)
library(MASS)
for(site.tag in 1:n.sites){
  for(visit.tag in 1:n.visit){
    C[visit.tag, site.tag] <- rcat(1, prop[1:n.species,site.tag])
  }
}


# Which species were present
#for(visit.tag in 1:n.visit){
#  for(site.tag in 1:n.sites){
 #   numb[visit.tag,site.tag] <-which(C[1:n.species,visit.tag,site.tag]==1)
 # }
#}


# Simulating the verified data
for(visit.tag in 1:n.visit){
  for(site.tag in 1:n.sites){
    Y[visit.tag,site.tag] <- rcat(1, prob = omega[C[visit.tag,site.tag], 1:(n.species+1)])
  }
}

#Sample only a handful of locations in the region
index <- sort(sample(1:n.sites,input$constants$sample.size, replace=FALSE))
#C_samp <- array(NA, dim=c(n.species, n.visit,input$constants$sample.size)) # Obs detect/non-detect
#Y_samp <- array(NA, dim=c(n.visit, input$constants$sample.size)) # Verified species with Machine Learning
C[,-index] <- NA
Y[,-index] <- NA
#for(site.tag in 1:length(index)){
 #   C_samp[,, site.tag] <- C[1:n.species,1:n.visit,index[site.tag]]
#}

#C_samp <- C[1:n.species,1:n.visit,index]
#for(site.tag in 1:length(index)){
#  Y_samp[, site.tag] <- Y[1:n.visit,index[site.tag]]
#}

  
locs_samp <- locs[index,]
cov_samp <- c(t(gridcov))[index]

# Which locations are selected
#for(visit.tag in 1:n.visit){
 # for(site.tag in 1:n.sites){
#    numb_with_NA[visit.tag,site.tag] <-  ifelse(sum(is.na(C[1:n.species,visit.tag,site.tag])) ==n.species, 0,1)
 # }
#}

#num_plot <- matrix(numb, nrow= dim[1], ncol=dim[2], byrow=F)
#num_with_NA_plot <- matrix(numb_with_NA, nrow= dim[1], ncol=dim[2], byrow=F)
if(input$plotdata$prop_plot ==TRUE){
data3 <- data.frame(x=locs$Var1,y=locs$Var2, maxprob = as.factor(max_prob),numb = as.factor(C[1,]), numb_na = as.factor(Y[1,]))
g1 <- ggplot(data3, aes(x,y))+
  geom_raster(aes(fill=numb))+
  scale_fill_manual(values = c("red", "yellow", "blue"), labels=c("S1","S2","S3"), guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("B) Verified species")

g2 <- ggplot(data3, aes(x,y))+
  geom_raster(aes(fill=maxprob))+
  scale_fill_manual(values = c("red", "yellow", "blue"), labels=c("S1","S2","S3"), guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("A) Species with max probability")

g3 <- ggplot(data3, aes(x,y))+
  geom_raster(aes(fill=numb_na))+
  scale_fill_manual(values = c("red", "yellow", "blue", "black"), labels=c("C1","C2","C3", "C4"), guide=guide_legend(reverse = TRUE))+
  theme_bw()+
  coord_fixed()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("C) Reported Species")

print(ggarrange(g2,g1,g3, nrow=2, ncol=2))
}

  
  dataout <- list(C=C,
                  #C_samp = C_samp,
                  Y=Y, 
                  #Y_samp = Y_samp,
                  locs=locs,
                  #locs_samp = locs_samp,
                  Lam=lambda,
                  omega=omega, 
                  cov=cov, 
                  #cov_samp=cov_samp,
                  prop=prop)
  return(dataout)
  
}

#genData(input)

nreplicates <- 50
sim <- pblapply(1:nreplicates, function(x){
  data <- genData(input, seed = x)
}, cl=1)


save(sim, file="simulated data1.RData")
