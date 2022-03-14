#This script runs an INLA within Nimble for the model in the main text
#Reader is referred to the vignette on multinomial logit models with INLA 
#https://rdrr.io/github/inbo/INLA/f/vignettes/multinomial.Rmd

#Load the data from "data formating.R"
load("yearly_data.RData")

#Packages needed to run the script
library(leaflet)
library(viridis)
library(sp)
library(rgdal)
library(tidyr)
library(raster)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
library(nimble)
library(INLA)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#Select the 2019 data
d <- sim[[1]] #sim is from loading yealy_data.RData
dim(unique(d[,c("latitude","longitude")]))

#transforming coordinates
sps <- SpatialPoints(d[,c("longitude","latitude")],
                     proj4string = CRS("+proj=utm +zone=28")
                     )
spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84"))

d[,c("long", "lat")] <- coordinates(spst)
head(d)

#Selecting the data that are in USA
init_plot <- SpatialPointsDataFrame(coords = cbind( d$longitude,d$latitude), data=d, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#plot(init_plot)
states <-raster::getData("GADM",country="United States",level=0)
points_in_norway <- raster::intersect(init_plot, states)
d <- points_in_norway@data[,-c(44,45)] #The data that are in USA

# Extract the altitude from Bioclim
r <- getData('worldclim', var='alt', res=10) 
pts <- SpatialPointsDataFrame(coords = d[,c("longitude","latitude")], data = d[,10:11], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
d$alt <- raster::extract(r, pts, df=T)

#Plot of data
ggplot(data = world)+
  geom_sf()+
  geom_point(data =d, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
  coord_sf(xlim = c(-135, -50), ylim = c(10.5, 58), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")+
  facet_wrap(~truth, ncol=2)+
  theme(legend.position = c(0.8,0.2))+
  labs(shape="Reported species")+
  labs(color="Reported species")

#Spread the data for multinomial regression, and remove data with errors
d$spread <- rep(1, nrow(d))
format_d <- d %>% spread(truth, spread, fill = 0, convert = FALSE)

indx_double <- vector("numeric", length=nrow(format_d))
for(site.tag in 1:nrow(format_d)){
indx_double[site.tag] <- sum(format_d[site.tag,44:46])
}
indx_delete <- which(indx_double==2)
d <- format_d[-indx_delete,] #data without formatting errors



#scale altitude
d$scale_alt <- scale(d$alt$alt,scale=FALSE)

#Dataframe with the structure for the multinomial regression in INLA
df = data.frame(cbind(d[,44:47], 
                      d$scale_alt, d$scale_alt,d$scale_alt,d$scale_alt,d$scale_alt,
                      d$num_identification_agreements,
                      d$num_identification_disagreements,
                      d$month,
                      d$longitude,
                      d$latitude,
                      as.numeric(as.factor(d$quality_grade))))
df <- as.matrix(df, rownames.force =TRUE)

# Fit the spatial multinomial regression
# Used to estimate the parameters of the true intensity
fit.inla <- function(df){
Data.structure = function(df){
  n=nrow(df)
  Data = matrix(NA, ncol = 12, nrow = n*3)
  for(i in 1:n){
    # simulated variable
    Data[((i-1)*3+1):(i*3), 1] = c(df[i,1], df[i,2], df[i,3])
    # alternative specific with altitude coeff
    Data[((i-1)*3+1):(i*3), 2] = c(df[i,4], df[i,4], df[i,4])
    # alternative specific with alternative coeff
    Data[((i-1)*3+1):(i*3), 3] = rep(df[i,10],3)
    Data[((i-1)*3+1):(i*3), 4] = rep(df[i,11],3)
    #Data[((i-1)*3+1):(i*3), 3:5] = diag(c(df$W.A[i], df$W.B[i], df$W.C[i]))
    # individual specific with alternative coeff
    Data[((i-1)*3+1):(i*3), 5] = rep(df[i,12],3)
    # choice situation index
    Data[((i-1)*3+1):(i*3), 6] = rep(i,3)
    # choice alternative index
    Data[((i-1)*3+1):(i*3), 7] = c(1, 2, 3)
    Data[((i-1)*3+1):(i*3), 8] = c(1, 2, 3)
    Data[((i-1)*3+1):(i*3), 9] = c(1, 2, 3)
    Data[((i-1)*3+1):(i*3), 10] = c(1, 2, 3)
    Data[((i-1)*3+1):(i*3), 11] = rep(df[i,15],3)
    Data[((i-1)*3+1):(i*3), 12] = c(1, 2, 3)
  }
  Data = data.frame(Data)
  names(Data) = c('Y', "alt","agree","disagree","month",'phi','alt.idx','alt.idx2','alt.idx3','alt.idx4','quality','alt.idx5' )
  return(Data)
}

#Prepare the data with the defined Data.structure function
new_df <- Data.structure(df)

#Mesh construction
coo <- cbind(rep(as.numeric(df[,13]),each=3), rep(as.numeric(df[,14]), each=3))
mesh <- INLA::inla.mesh.2d(loc=coo,
                     max.edge=c(2,5),
                     cutoff = 0.01
)
#Plot mesh
#plot(mesh)
#points(coo, col="red")

#Building the SPDE model on the mesh
spde=INLA::inla.spde2.matern(mesh=mesh,
                       alpha=2,
                       constr=TRUE)

#index set
indexs <- INLA::inla.spde.make.index("s", spde$n.spde)
#lengths(indexs)

#Projection matrix
A <- INLA::inla.spde.make.A(mesh=mesh, loc=coo)
  
#stack data
stk.e <- INLA::inla.stack(
  tag="est",
  data=list(y=as.numeric(new_df$Y)),
  A=list(1,A),
  effects=list(data.frame(b0=1, 
                          altitude=as.numeric(new_df$alt),
                          agree=as.numeric(new_df$agree),
                          disagree = as.numeric(new_df$disagree),
                          month=as.factor(new_df$month),
                          phi=as.numeric(new_df$phi),
                          quality=as.factor(new_df$quality),
                          alt.idx = as.numeric(new_df$alt.idx),
                          alt.idx2=as.numeric(new_df$alt.idx2),
                          alt.idx3=as.numeric(new_df$alt.idx3),
                          alt.idx4=as.numeric(new_df$alt.idx4),
                          alt.idx5=as.numeric(new_df$alt.idx5)),
                          s=indexs)
)

#INLA formula
formula <- y ~ 0+b0+
  f(phi, initial = -10, fixed=TRUE)+
  f(alt.idx, agree, fixed = T, constr = T)+
  f(alt.idx2, disagree, fixed = T, constr = T)+
  f(alt.idx3, factor(month), fixed = T, constr = T)+
  f(alt.idx4, altitude, fixed = T, constr = T)+
  f(alt.idx5, factor(quality), fixed = T, constr = T)+
  f(s, model=spde)
  
  #Fit the model
res <- INLA::inla(formula,
            family="poisson",
            control.family = list(link="log"),
            data=inla.stack.data(stk.e),
            control.predictor = list(
              compute=TRUE, link=1,
              A=inla.stack.A(stk.e)
            ),
            control.inla = list(int.strategy="eb"))

  #selecting the fitted values of lambda
idx <- INLA::inla.stack.index(stk.e,'est')$data
fitted_val <- matrix(res$summary.fitted.values[idx,"mean"], ncol=3, byrow=FALSE)
return(fitted_val)
}
#fit.inla(df)

# Making the INLA work in NIMBLE
nimbleINLA <- nimbleRcall(
  prototype = function(
    df=double(2)
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)
#Compile nimble model
CnimbleINLA <- compileNimble(nimbleINLA)
CnimbleINLA(df)

#Code to estimate classification probabilities
#The lambda's are estimated from INLA and 
#used to estimate classification probabilities
code <- nimbleCode({

#Linear predictor
  lambda[1:n.sites,1:n.species] <- nimbleINLA(df[1:n.sites, 1:15])
  
  # Proportion for the multinomial distribution
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      prop[site.tag,spe.tag] <- lambda[site.tag, spe.tag]/sum(lambda[site.tag, 1:n.species])
    }
  }
  
  #Prior for the classification probabilities
  for(spe.tag in 1:n.species){
    for(i in 1:n.reported){
      alpha[spe.tag, i] ~ dexp(1)
    }
  }
  
  # Confusion Matrix for the Species
  for(spe.tag in 1:n.species){
    omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
  }
  
  # Verified observations 
  for(site.tag in 1:n.sites){
    C[site.tag] ~ dcat(prop[site.tag,1:n.species])
  }
  
  # Verified observations from ML classifications
  for(site.tag in 1:n.sites){
    Y[site.tag] ~ dcat(omega[C[site.tag],1:n.reported])
  }
})

# Data for NIMBLE
#verified data
C <- vector("numeric", length=nrow(d))
for(i in 1:nrow(d)){
  C[i] <- which(c(d[i,44:46])==1)
}

# Reported observations
Y <- as.numeric(as.factor(d$taxon_species_name))

#Setting parameters
data <- list(C,Y,df)
n.species <- length(unique(c(data[[1]])))
n.reported <- length(unique(c(data[[2]])))
n.sites <- length(data[[1]])

# Constants for the MCMC model
const <- list(n.sites=n.sites,
              n.species=n.species,
              n.reported=n.reported)

# Data for the Model
idm_data <- list(C = data[[1]], 
                 Y = data[[2]],
                 df=data[[3]])

# Initials for the model
alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.reported))
for(spe.tag in 1:const$n.species){
  for(k in 1:(const$n.reported)){
    alpha[spe.tag,k] <- 1
  }
}
omega <- matrix(NA, nrow=const$n.species,
                ncol = (const$n.reported)) 
for(spe.tag in 1:(const$n.species)){
  omega[spe.tag,]<- MCMCpack::rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.reported)])
}

idm_inits <- function(){list(omega=omega,
                             alpha = alpha)
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
                    inits = initsList)

# Create the model in C
Cmwtc <- compileNimble(mwtc, 
                       showCompilerOutput = FALSE)

mcmcconf <- configureMCMC(Cmwtc, monitors = c("omega"), 
                          enableWAIC =TRUE)
Rmcmc <- buildMCMC(mcmcconf)

# Compile 
cmcmc <- compileNimble(Rmcmc, 
                       project = Cmwtc,
                       resetFunctions = TRUE)

# Run the MCMC
mcmc.out <- runMCMC(cmcmc, 
                    niter = 5000,
                    nchains = 3,
                    nburnin = 2000,
                    inits = initsList,
                    setSeed = TRUE, 
                    samples=TRUE, 
                    samplesAsCodaMCMC = TRUE, 
                    summary = TRUE, 
                    WAIC = FALSE)
save(mcmc.out, file="mcmc_result.RData")

#Check convergence of results
#library(ggmcmc)
#mcmclist <- ggs(mcmc.out$samples)
#ggs_density(mcmclist, family = "omega")
#ggs_density(mcmclist, family = "beta")
#ggs_traceplot(mcmclist, family = "beta")
#ggs_traceplot(mcmclist, family = "omega")

output <- mcmc.out$summary$all.chains
omega <- matrix(output[,2], 
                nrow=const$n.species, 
                ncol=(const$n.reported), 
                byrow=FALSE)

#save classification probabilities
save(omega, file="omega.RData")

#################################
# Prediction for 2019
#################################
#Fitted values for the true intensity of verified species
lambda <- fit.inla(df)

#select the data that needs ID in 2019
indx_nid <- which(d$quality_grade=="needs_id")

#Function to estimate the percentage of correct identification 
#before and after prediction
proportions_estimator <- function(d, lambda,indx_nid){
lambda_nid <- lambda[indx_nid,] #lambda for the selected subset that needs ID 

prop_nid <- t(apply(lambda_nid,1,proportions)) #estimate proportions

 #extract reported species  
reported_species <- as.numeric(as.factor(d[indx_nid,]$taxon_species_name))

#Estimating posterior probabilities 
post_prob <- prop_true_sp <- matrix(NA, ncol =3, nrow=nrow(prop_nid))
for(site.tag in 1:nrow(prop_nid)){
  for(spe.tag in 1:3){
  post_prob[site.tag,spe.tag] <- omega[spe.tag,reported_species[site.tag]]*prop_nid[site.tag, spe.tag]
}
for(site.tag in 1:nrow(prop_nid)){
  for(spe.tag in 1:3){
    prop_true_sp[site.tag,spe.tag] <- post_prob[site.tag,spe.tag]/(rowSums(post_prob)[site.tag])
  }
}

#species with highest posterior probabilities
pred_true_sps <- vector("numeric", length=nrow(prop_nid))
for(site.tag in 1:nrow(prop_nid)){
  pred_true_sps[site.tag] <- which.max(prop_true_sp[site.tag,])
}

  #Extracting verified data
C <- vector("numeric", nrow(prop_nid))
for(i in 1:nrow(prop_nid)){
  C[i] <- which(c(d[indx_nid[i],44:46])==1)
}

  #estimating percentages
percentage_correct_before_pred <- mean(C==reported_species)
percentage_correct_after_pred <- mean(pred_true_sps==C)

ret <- data.frame(reported_species=reported_species,
                  pred_true_sps=pred_true_sps,
                  true_species = C)
return(list(ret=ret,
            percentage_correct_before_pred=percentage_correct_before_pred,
              percentage_correct_after_pred=percentage_correct_after_pred,
            prop_true_sp=prop_true_sp))
}

  #estimating the needed percentages for the data that need ID in 2019 and save result
pred_true_species_2019_needs_19 <- proportions_estimator(d, lambda,indx_nid)
save(pred_true_species_2019_needs_19, file="pred_needs_19.RData")
  
####################################
#Prediction for 2020
####################################
  
d1 <- sim[[2]] #select the 2020 data
dim(unique(d1[,c("latitude","longitude")]))

#transforming coordinates
sps1 <- SpatialPoints(d1[,c("longitude","latitude")],
                     proj4string = CRS("+proj=utm +zone=28")
)
spst1 <- spTransform(sps1, CRS("+proj=longlat +datum=WGS84"))

d1[,c("long", "lat")] <- coordinates(spst1)
head(d1)

  #Extracting the data for USA
init_plot1 <- SpatialPointsDataFrame(coords = cbind( d1$longitude,d1$latitude), data=d1, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
states <-raster::getData("GADM",country="United States",level=1)
selected_stated <- states[states$NAME_1=="Florida"|states$NAME_1=="Montana"| states$NAME_1=="Minnesota"| states$NAME_1=="Wisconsin"| states$NAME_1=="Illinois"| states$NAME_1=="Indiana",]
points_in_norway1 <- raster::intersect(init_plot1, selected_stated)
d1 <- points_in_norway1@data[,-c(47,53)] #The data in 2020 from the selected states
 
  #extract the altitude data
r <- getData('worldclim', var='alt', res=10) 
pts1 <- SpatialPointsDataFrame(coords = d1[,c("longitude","latitude")], data = d1[,10:11], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
d1$alt <- raster::extract(r, pts1, df=T)

#remove other danauus, since the 2019 data does not have that category
#for easy prediction for the 2020
d1 <- d1 %>% 
  dplyr::filter(taxon_species_name != "other danaus")%>%
  dplyr::filter(!is.na(alt$alt))

#Plot of data
ggplot(data = world)+
  geom_sf()+
  geom_point(data =d1, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
  coord_sf(xlim = c(-135, -50), ylim = c(10.5, 58), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")+
  facet_wrap(~truth, ncol=2)+
  theme(legend.position = c(0.8,0.2))+
  labs(shape="Reported species")+
  labs(color="Reported species")

#Spread the data for multinomial regression
d1$spread1 <- rep(1, nrow(d1))
format_d1 <- d1 %>% 
  tidyr::spread(truth, spread1, fill = 0)#, convert = TRUE)
indx_double1 <- vector("numeric", length=nrow(format_d1))
for(site.tag in 1:nrow(format_d1)){
  indx_double1[site.tag] <- sum(format_d1[site.tag,52:54])
}
indx_delete1 <- which(indx_double1!=1)
format_d1 <- format_d1[-indx_delete1,]

#scale altitude
d1 <- format_d1
d1$scale_alt <- d1$alt$alt - mean(d$alt$alt, na.rm=T)

  
#Data structure for multinomial
df1 = data.frame(cbind(d1[,52:54], 
                      d1$scale_alt, d1$scale_alt,d1$scale_alt,d1$scale_alt,d1$scale_alt,
                      d1$num_identification_agreements,
                      d1$num_identification_disagreements,
                      d1$month,
                      d1$longitude,
                      d1$latitude,
                      as.numeric(as.factor(d1$quality_grade))))
df_pred <- as.matrix(df1, rownames.force =TRUE)

#Fit INLA with 2019 data for estimation and the 2020 data for prediction. 
fit.inla_prediction <- function(df, df_pred){
  Data.structure = function(df){
    n=nrow(df)
    Data = matrix(NA, ncol = 12, nrow = n*3)
    for(i in 1:n){
      # simulated variable
      Data[((i-1)*3+1):(i*3), 1] = c(df[i,1], df[i,2], df[i,3])
      # alternative specific with generic coeff
      Data[((i-1)*3+1):(i*3), 2] = c(df[i,4], df[i,4], df[i,4])
      # alternative specific with alternative coeff
      Data[((i-1)*3+1):(i*3), 3] = rep(df[i,10],3)
      Data[((i-1)*3+1):(i*3), 4] = rep(df[i,11],3)
      #Data[((i-1)*3+1):(i*3), 3:5] = diag(c(df$W.A[i], df$W.B[i], df$W.C[i]))
      # individual specific with alternative coeff
      Data[((i-1)*3+1):(i*3), 5] = rep(df[i,12],3)
      # choice situation index
      Data[((i-1)*3+1):(i*3), 6] = rep(i,3)
      # choice alternative index
      Data[((i-1)*3+1):(i*3), 7] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 8] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 9] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 10] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 11] = rep(df[i,15],3)
      Data[((i-1)*3+1):(i*3), 12] = c(1, 2, 3)
    }
    Data = data.frame(Data)
    names(Data) = c('Y', "alt","agree","disagree","month",'phi','alt.idx','alt.idx2','alt.idx3','alt.idx4','quality','alt.idx5' )
    return(Data)
  }
  
  Data.structure1 = function(df1){
    n=nrow(df1)
    Data = matrix(NA, ncol = 12, nrow = n*3)
    for(i in 1:n){
      # simulated variable
      Data[((i-1)*3+1):(i*3), 1] = c(df1[i,1], df1[i,2], df1[i,3])
      # alternative specific with generic coeff
      Data[((i-1)*3+1):(i*3), 2] = c(df1[i,4], df1[i,4], df1[i,4])
      # alternative specific with alternative coeff
      Data[((i-1)*3+1):(i*3), 3] = rep(df1[i,9],3)
      Data[((i-1)*3+1):(i*3), 4] = rep(df1[i,10],3)
      #Data[((i-1)*3+1):(i*3), 3:5] = diag(c(df$W.A[i], df$W.B[i], df$W.C[i]))
      # individual specific with alternative coeff
      Data[((i-1)*3+1):(i*3), 5] = rep(df1[i,11],3)
      # choice situation index
      Data[((i-1)*3+1):(i*3), 6] = rep(i,3)
      # choice alternative index
      Data[((i-1)*3+1):(i*3), 7] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 8] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 9] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 10] = c(1, 2, 3)
      Data[((i-1)*3+1):(i*3), 11] = rep(df1[i,14],3)
      Data[((i-1)*3+1):(i*3), 12] = c(1, 2, 3)
    }
    Data = data.frame(Data)
    names(Data) = c('Y', "alt","agree","disagree","month",'phi','alt.idx','alt.idx2','alt.idx3','alt.idx4','quality','alt.idx5' )
    return(Data)
  }
  #round(head(Data.structure(df)),3)
  new_df <- Data.structure(df)
  new_df1 <- Data.structure1(df_pred)
  
  #Mesh construction
  coo <- cbind(rep(as.numeric(df[,13]),each=3), rep(as.numeric(df[,14]), each=3))
  mesh <- INLA::inla.mesh.2d(loc=coo,
                             max.edge=c(2,5),
                             cutoff = 0.01
  )

 spde = INLA::inla.spde2.matern(mesh=mesh,
                          alpha=2,
                          constr=TRUE)
  
  #index set
  indexs <- INLA::inla.spde.make.index("s", spde$n.spde)
  #lengths(indexs)
  
  #Projection matrix
  A <- INLA::inla.spde.make.A(mesh=mesh, loc=coo)
  #stack data
  stk.e <- INLA::inla.stack(
    tag="est",
    data=list(y=as.numeric(new_df$Y)),
    A=list(1,A),
    effects=list(data.frame(b0=1, 
                            altitude=as.numeric(new_df$alt),
                            agree=as.numeric(new_df$agree),
                            disagree = as.numeric(new_df$disagree),
                            month=as.factor(new_df$month),
                            phi=as.numeric(new_df$phi),
                            quality=as.factor(new_df$quality),
                            alt.idx = as.numeric(new_df$alt.idx),
                            alt.idx2=as.numeric(new_df$alt.idx2),
                            alt.idx3=as.numeric(new_df$alt.idx3),
                            alt.idx4=as.numeric(new_df$alt.idx4),
                            alt.idx5=as.numeric(new_df$alt.idx5)),
                 s=indexs)
  )
  
  
  #Prediction points
  pts.pred <- cbind(rep(as.numeric(df_pred[,13]),each=3), rep(as.numeric(df_pred[,14]), each=3))
  A.pred <- inla.spde.make.A(mesh=mesh, loc=pts.pred)
  
  stk.p <- INLA::inla.stack(
    tag="pred",
    data=list(y=NA),
    A=list(1,A.pred),
    effects=list(data.frame(b0=1, 
                            altitude=as.numeric(new_df1$alt),
                            agree=as.numeric(new_df1$agree),
                            disagree = as.numeric(new_df1$disagree),
                            month=as.factor(new_df1$month),
                            phi=as.numeric(new_df1$phi),
                            quality=as.factor(new_df1$quality),
                            alt.idx = as.numeric(new_df1$alt.idx),
                            alt.idx2=as.numeric(new_df1$alt.idx2),
                            alt.idx3=as.numeric(new_df1$alt.idx3),
                            alt.idx4=as.numeric(new_df1$alt.idx4),
                            alt.idx5=as.numeric(new_df1$alt.idx5)),
                 s=indexs)
  )
  
  joint.stk <- inla.stack(stk.e, stk.p)
  
  #INLA formula
  formula <- y ~ 0+b0+
    f(phi, initial = -10, fixed=TRUE)+
    f(alt.idx, agree, fixed = T, constr = T)+
    f(alt.idx2, disagree, fixed = T, constr = T)+
    f(alt.idx3, factor(month), fixed = T, constr = T)+
    f(alt.idx4, altitude, fixed = T, constr = T)+
    f(alt.idx5, factor(quality), fixed = T, constr = T)+
    f(s, model=spde, group = s.group)
  
  #Fitting the model and saving results
  res <- INLA::inla(formula,
                    family="poisson",
                    control.family = list(link="log"),
                    data=inla.stack.data(joint.stk, spde=spde),
                    control.predictor = list(
                      compute=TRUE, link=1,
                      A=inla.stack.A(joint.stk)
                    ),
                    control.compute = list(config=TRUE),
                    control.inla = list(int.strategy="eb"))
  save(res, file="prediction.RData")

  #select the fitted values for the prediction
  idx <- INLA::inla.stack.index(stk.p,'pred')$data
  fitted_val <- matrix(res$summary.fitted.values[idx,"mean"], ncol=3, byrow=FALSE)
  return(fitted_val)
}

  #Predicted fitted values for 2020
lambda_pred <- fit.inla_prediction(df,df_pred)

  ## estimating percentage of correct classification
proportions_estimator_pred <- function(d1, lambda_pred,indx_nid_pred){
  lambda_nid <- lambda_pred[indx_nid_pred,] #indexing the lambda 
  prop_nid <- t(apply(lambda_nid,1,proportions)) #proportions
 
  #extracting reported species
  reported_species <- as.numeric(as.factor(d1[indx_nid_pred,]$taxon_species_name))
  
    #estimating posterior probability
  post_prob <- prop_true_sp <- matrix(NA, ncol =3, nrow=nrow(prop_nid))
  for(site.tag in 1:nrow(prop_nid)){
    for(spe.tag in 1:3){
      post_prob[site.tag,spe.tag] <- omega[spe.tag,reported_species[site.tag]]*prop_nid[site.tag, spe.tag]
    }
  }
  for(site.tag in 1:nrow(prop_nid)){
    for(spe.tag in 1:3){
      prop_true_sp[site.tag,spe.tag] <- post_prob[site.tag,spe.tag]/(rowSums(post_prob)[site.tag])
    }
  }
  
  #predicted true species
  pred_true_sps <- vector("numeric", length=nrow(prop_nid))
  for(site.tag in 1:nrow(prop_nid)){
    pred_true_sps[site.tag] <- which.max(prop_true_sp[site.tag,])
  }
  
  # Verified species
  C <- vector("numeric", nrow(prop_nid))
  for(i in 1:nrow(prop_nid)){
    C[i] <- which(d1[indx_nid_pred[i],52:54]==1)
  }
  
  percentage_correct_before_pred <- mean(C==reported_species)
  percentage_correct_after_pred <- mean(pred_true_sps==C)
  
  ret <- data.frame(reported_species=reported_species,
                    pred_true_sps=pred_true_sps,
                    true_species = C)
  return(list(ret=ret,
              percentage_correct_before_pred=percentage_correct_before_pred,
              percentage_correct_after_pred=percentage_correct_after_pred,
              prop_true_sp=prop_true_sp))
}

  #2020 data that need ID
indx_nid_pred <- which(d1$quality_grade=="needs_id")
pred_true_species_2020_needs_20 <- proportions_estimator_pred(d1, lambda_pred,indx_nid_pred)
save(pred_true_species_2020_needs_20, file="needs_20.RData")

  #2020 data with research grade
indx_nid_pred <- which(d1$quality_grade=="research")
pred_true_species_2020_research_20 <- proportions_estimator_pred(d1, lambda_pred,indx_nid_pred)
save(pred_true_species_2020_research_20, file="research_20.RData")




