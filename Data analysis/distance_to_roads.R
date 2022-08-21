#Function that estimates distance to road
dist_to_road <- function(data){
  
  #packages
  library(tigris)
  library(dplyr)
  library(sf)
  library(sp)
  library(raster)
  library(stringr)#install.packages("usmap")
  library(usmap)
  
groups <- data %>%
  arrange(NAME_1, NAME_2)

#data[viceroy$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"
groups[grep("^Saint", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Saint", groups$NAME_2),]$NAME_2, "Saint", "St.") 
groups[grep("^St.e Genevieve", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^St.e Genevieve", groups$NAME_2),]$NAME_2, "^St.e Genevieve", "Ste. Genevieve") 
groups[grep("^Lake Hurron", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Hurron", groups$NAME_2),]$NAME_2, "^Lake Hurron", "Huron") 
groups[grep("^Lake St. Clair", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake St. Clair", groups$NAME_2),]$NAME_2, "^Lake St. Clair", "St. Clair") 
groups[grep("^Lake Superior", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Superior", groups$NAME_2),]$NAME_2, "^Lake Superior", "Lake") 
groups[groups$NAME_2=="De Kalb",]$NAME_2 <- "DeKalb"
groups[groups$NAME_2=="La Salle",]$NAME_2 <- "LaSalle"
groups[groups$NAME_2=="Debaca",]$NAME_2 <- "De Baca"
groups[groups$NAME_2=="Mc Kean",]$NAME_2 <- "McKean"
groups[groups$NAME_2=="Alexandria",]$NAME_2 <- "Alexandria city"
groups[groups$NAME_2=="Charlottesville",]$NAME_2 <- "Charlottesville city"
groups[groups$NAME_2=="Chesapeake",]$NAME_2 <- "Accomack"
groups[groups$NAME_2=="Falls Church",]$NAME_2 <- "Falls Church city"
groups[groups$NAME_2=="Fredericksburg",]$NAME_2 <- "Fredericksburg city"
groups[groups$NAME_2=="Hampton",]$NAME_2 <- "Hampton city"
groups[groups$NAME_2=="Harrisonburg",]$NAME_2 <- "Harrisonburg city"
groups[groups$NAME_2=="Lynchburg",]$NAME_2 <- "Lynchburg city"
groups[groups$NAME_2=="Manassas",]$NAME_2 <- "Manassas city"
groups[groups$NAME_2=="Newport News",]$NAME_2 <- "Newport News city"
groups[groups$NAME_2=="Norton",]$NAME_2 <- "Norton city"
groups[groups$NAME_2=="Poquoson",]$NAME_2 <- "Poquoson city"
groups[groups$NAME_2=="Portsmouth",]$NAME_2 <- "Portsmouth city"
groups[groups$NAME_2=="Radford",]$NAME_2 <- "Radford city"
groups[groups$NAME_2=="Richmond",]$NAME_2 <- "Richmond city"
groups[groups$NAME_2=="Virginia Beach",]$NAME_2 <- "Virginia Beach city"
groups[groups$NAME_2=="Williamsburg",]$NAME_2 <- "Williamsburg city"
groups[groups$NAME_2=="Winchester",]$NAME_2 <- "Winchester city"
groups[groups$NAME_2=="Waynesboro",]$NAME_2 <- "Waynesboro city"
groups[groups$NAME_2=="Richmond city",]$NAME_2 <- "Richmond"
groups[groups$NAME_2=="Fairbanks North Star",]$NAME_2 <- "Fairbanks North Star Borough"
groups[grep("^Lake Michigan", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Michigan", groups$NAME_2),]$NAME_2, "Lake Michigan", "Lake") #groups[groups$NAME_2=="Saint Clair",]$NAME_2 <- "St. Clair"
groups[grep("^Lake Ontario", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Ontario", groups$NAME_2),]$NAME_2, "Lake Ontario", "Ontario")
groups[grep("^Lake Erie", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Erie", groups$NAME_2),]$NAME_2, "Lake Erie", "Erie")

#extracting names of states
states_and_county <- groups%>%
  dplyr::select("NAME_1")
names_states <- unique(states_and_county$NAME_1)

x <- list()
for(states in 1:length(names_states)){
   print(states)
  y <- list()
  subset_data <- groups%>%
    filter(NAME_1 == names_states[states])

  #remove Lake from Wisconsin
  if(names_states[states]=="Wisconsin"){
    subset_data <- subset_data%>%
      filter(NAME_2 !="Lake") 
  }else{
    subset_data <- subset_data 
  }
  
  if(names_states[states]=="Massachusetts"){
    subset_data[subset_data$NAME_2=="Norfolk",]$NAME_2 <- "Norfolk"
    subset_data[subset_data$NAME_2=="Suffolk",]$NAME_2 <- "Suffolk"
  }else{
    if(names_states[states]=="Virginia"){
    subset_data[subset_data$NAME_2=="Norfolk",]$NAME_2 <- "Norfolk city"
    subset_data[subset_data$NAME_2=="Suffolk",]$NAME_2 <- "Suffolk city"
    }else{
      subset_data <- subset_data
    }
  }
  
  names_county <- unique(subset_data$NAME_2)
  
  #Add parish to county names
  if(names_states[states]=="Louisiana"){
    new_names_county <- paste(names_county, "parish", sep = " ") 
  }else{
    new_names_county <- names_county
  }
  
  
  fips_codes <- vector("character", length(names_county))
  
  for(i in 1:length(names_county)){
  fips_codes[i] = fips(state = names_states[states],
                    county = new_names_county[i])%>%
    stringr::str_sub(.,3,5)
  }
  if(length(names_county) != length(fips_codes)) stop("Check the county name and FIPS codes")
  
  for(county in 1:length(names_county)){
    print(county)
    if(names_county[county]=="Ada"){
      roads1 <- roads(names_states[states],county= "001")}else{
    roads1 <- roads(names_states[states],county= fips_codes[county])
      }
    subset_data1 <- subset_data%>%
      filter(NAME_2 == names_county[county])
    coordinates <- cbind(subset_data1$longitude, subset_data1$latitude)
    if(nrow(coordinates > 0)){
    dist_grid <- st_as_sf(SpatialPointsDataFrame(coords = coordinates ,data=subset_data1,proj4string = crs(roads1)))
    dists_0 <- st_distance(dist_grid,roads1)
    dist2road <- apply(dists_0,1,min)
    }else{
      dist2road = 0 
    }
    subset_data1$dist_road <- dist2road
    y[[county]] <- subset_data1
  }
  x[[states]] <- y
}
all_data <- purrr::flatten(x)
returned_dataframe <- do.call("rbind", all_data)
return(returned_dataframe)
}

dist_to_road1 <- function(data){
  
  #packages
  library(tigris)
  library(dplyr)
  library(sf)
  library(sp)
  library(raster)
  library(stringr)#install.packages("usmap")
  library(usmap)
  
  groups <- data %>%
    arrange(NAME_1, NAME_2)
  
  #data[viceroy$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"
  groups[grep("^Saint", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Saint", groups$NAME_2),]$NAME_2, "Saint", "St.") 
  groups[grep("^St.e Genevieve", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^St.e Genevieve", groups$NAME_2),]$NAME_2, "^St.e Genevieve", "Ste. Genevieve") 
  groups[grep("^Lake Hurron", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Hurron", groups$NAME_2),]$NAME_2, "^Lake Hurron", "Huron") 
  groups[grep("^Lake St. Clair", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake St. Clair", groups$NAME_2),]$NAME_2, "^Lake St. Clair", "St. Clair") 
  groups[grep("^Lake Superior", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Superior", groups$NAME_2),]$NAME_2, "^Lake Superior", "Lake") 
  groups[grep("^Lake Michigan", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Michigan", groups$NAME_2),]$NAME_2, "Lake Michigan", "Lake") #groups[groups$NAME_2=="Saint Clair",]$NAME_2 <- "St. Clair"
  groups[grep("^Lake Ontario", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Ontario", groups$NAME_2),]$NAME_2, "Lake Ontario", "Ontario")
  groups[grep("^Lake Erie", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Lake Erie", groups$NAME_2),]$NAME_2, "Lake Erie", "Erie")
  groups[grep("^Debaca", groups$NAME_2),]$NAME_2 <- str_replace_all(groups[grep("^Debaca", groups$NAME_2),]$NAME_2, "Debaca", "De baca")
  groups[groups$NAME_2=="De Kalb",]$NAME_2 <- "DeKalb"
  #groups[groups$NAME_2=="Shannon",]$NAME_2 <- "Bennett"
  #extracting names of states
  states_and_county <- groups%>%
    dplyr::select("NAME_1")
  names_states <- unique(states_and_county$NAME_1)
  
  #print
  
  x <- list()
  for(states in 1:length(names_states)){
    print(states)
    y <- list()
    subset_data <- groups%>%
      filter(NAME_1 == names_states[states])
    
    #remove Lake from Wisconsin
    if(names_states[states]=="Wisconsin"){
      subset_data <- subset_data%>%
        filter(NAME_2 !="Lake") 
    }else{
      subset_data <- subset_data 
    }
    
    #if(names_states[states]=="Massachusetts"){
      #subset_data[subset_data$NAME_2=="Norfolk",]$NAME_2 <- "Norfolk"
      #subset_data[subset_data$NAME_2=="Suffolk",]$NAME_2 <- "Suffolk"
    #}
    #else{
     # if(names_states[states]=="Virginia"){
     #   subset_data[subset_data$NAME_2=="Norfolk",]$NAME_2 <- "Norfolk city"
     #   subset_data[subset_data$NAME_2=="Suffolk",]$NAME_2 <- "Suffolk city"
     # }else{
     #   subset_data <- subset_data
     # }
    #}
    
    names_county <- unique(subset_data$NAME_2)
    
    #Add parish to county names
    if(names_states[states]=="Louisiana"){
      new_names_county <- paste(names_county, "parish", sep = " ") 
    }else{
      new_names_county <- names_county
    }
    
    
    fips_codes <- vector("character", length(names_county))
    
    for(i in 1:length(names_county)){
      fips_codes[i] = fips(state = names_states[states],
                           county = new_names_county[i])%>%
        stringr::str_sub(.,3,5)
    }
    if(length(names_county) != length(fips_codes)) stop("Check the county name and FIPS codes")
    
    for(county in 1:length(names_county)){
      print(county)
      if(names_county[county]=="Ada"){
        roads1 <- roads(names_states[states],county= "001")}else{
          roads1 <- roads(names_states[states],county= fips_codes[county])
        }
      subset_data1 <- subset_data%>%
        filter(NAME_2 == names_county[county])
      coordinates <- cbind(subset_data1$longitude, subset_data1$latitude)
      if(nrow(coordinates > 0)){
        dist_grid <- st_as_sf(SpatialPointsDataFrame(coords = coordinates ,data=subset_data1,proj4string = crs(roads1)))
        dists_0 <- st_distance(dist_grid,roads1)
        dist2road <- apply(dists_0,1,min)
      }else{
        dist2road = 0 
      }
      subset_data1$dist_road <- dist2road
      y[[county]] <- subset_data1
    }
    x[[states]] <- y
  }
  all_data <- purrr::flatten(x)
  returned_dataframe <- do.call("rbind", all_data)
  return(returned_dataframe)
}

#distance_to_road_data <- dist_to_road(d)

#fips_info("51740")

#All points



load("yearly_data_new.RData")
#source("distance_to_roads.R")

#training_data_with_distance <- dist_to_road1(sim[[1]])

all_data_with_distance <- rbind(distance_to_road[[1]], distance_to_road[[2]])%>%
  dplyr::filter(NAME_1 %in% c(  "Texas"))%>%
  dplyr::group_by(user_id)%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(count = n(),
                previous_obs = row_number() - 1)%>%
  #dplyr::arrange(desc(count))%>%
  #dplyr::filter(count >= 300 & count <= 900)%>%
  dplyr::ungroup()%>%
  dplyr::mutate(experience  = ifelse(count < 50, 1, 
                                     ifelse(count >50 & count < 100, 2, 
                                            ifelse(count > 100 & count < 200, 3,
                                                   ifelse(count > 300 & count < 400, 4, 
                                                          ifelse(count > 400 & count < 500, 5,  6))))))





# training_data_with_distance <- distance_to_road[[1]]%>%
#   dplyr::filter(NAME_1 %in% c(  "Texas"))%>%
#   dplyr::group_by(user_id)%>%
#   dplyr::arrange(id)%>%
#   dplyr::mutate(count = n(),
#                 previous_obs = row_number() - 1)%>%
#   #dplyr::arrange(desc(count))%>%
#   #dplyr::filter(count >= 300 & count <= 900)%>%
#   dplyr::ungroup()%>%
#   dplyr::mutate(experience  = ifelse(count < 100, 1, 
#                                      ifelse(count >100 & count < 200, 2, 
#                                             ifelse(count > 300 & count < 400, 3,
#                                                    ifelse(count > 400 & count < 500, 4, 
#                                                           ifelse(count > 500 & count < 600, 5,  6))))))
#   #dplyr::filter(year %in% c("2018", "2019"))%>%
# 
# table(training_data_with_distance$truth, training_data_with_distance$taxon_species_name)
# #validation_data_with_distance <- dist_to_road1(sim[[2]])
# validation_data_with_distance <- distance_to_road[[2]]%>%
#   dplyr::filter(NAME_1 %in% c(  "Texas"))%>%
#   dplyr::group_by(user_id)%>%
#   dplyr::arrange(id)%>%
#   dplyr::mutate(count = n(),
#                 previous_obs = row_number() - 1)%>%
#   #dplyr::arrange(desc(count))%>%
#   #dplyr::filter(count >= 300 & count <= 900)%>%
#   dplyr::ungroup()%>%
#   dplyr::mutate(experience  = ifelse(count < 100, 1, 
#                                      ifelse(count >100 & count < 200, 2, 
#                                             ifelse(count > 300 & count < 400, 3,
#                                                    ifelse(count > 400 & count < 500, 4, 
#                                                           ifelse(count > 500 & count < 600, 5, 6))))))
# #"Florida", "Alabama", "Georgia", 
# table(validation_data_with_distance$truth, validation_data_with_distance$taxon_species_name)
#training_data_with_distance <- training_data_with_distance[1:400,]
#validation_data_with_distance <- validation_data_with_distance[1:200,]
#prediction_data_with_distance <- dist_to_road1(sim[[3]])

validation_data_with_distance <- all_data_with_distance%>%
  filter(year == 2020)

training_data_with_distance <- all_data_with_distance%>%
  filter(year!= 2020)

Y <- c(training_data_with_distance$taxon_species_name, validation_data_with_distance$taxon_species_name)
C <- c(training_data_with_distance$truth, validation_data_with_distance$truth)

index_for_splitting <- seq(nrow(training_data_with_distance)+1, 
                           nrow(training_data_with_distance)+nrow(validation_data_with_distance))

Validation_C <- as.numeric(as.factor(validation_data_with_distance$truth))
Validation_Y <- as.numeric(as.factor(validation_data_with_distance$taxon_species_name))

validation_indices_for_mismatch = index_for_splitting[which(Validation_C != Validation_Y)]
validation_correct_match = index_for_splitting[which(Validation_C == Validation_Y)]
validation_mismatchC = C[validation_indices_for_mismatch]
validation_mismatchY = Y[validation_indices_for_mismatch]
validation_correctmatchC = C[validation_correct_match]
validation_correctmatchY = Y[validation_correct_match]
#table(Validation_C,Validation_Y)
C[index_for_splitting] <- NA
#Y[1:n.visit, index_for_splitting] <- NA

#covariates
altitude <- c(training_data_with_distance$alt, validation_data_with_distance$alt)/1000
distance <- c(training_data_with_distance$dist_road, validation_data_with_distance$dist_road)/1000
precipitation <- c(training_data_with_distance$bio12, validation_data_with_distance$bio12)/1000
mean_temperature <- c(training_data_with_distance$bio1, validation_data_with_distance$bio1)/10

#unique user
users <- as.numeric(as.factor( c(training_data_with_distance$user_id, validation_data_with_distance$user_id)))

#experience
#experience <- c(training_data_with_distance$previous_obs, validation_data_with_distance$previous_obs)
experience <- c(training_data_with_distance$experience, validation_data_with_distance$experience)
no_previous_obs <- c(training_data_with_distance$previous_obs, validation_data_with_distance$previous_obs)


dataout <- list(C=C, #Verified data
                Y=Y, # Reported data
                altitude=altitude, #Mean intensity
                distance=distance, # Confusion matrix
                users=users, # covariate
                #cov_omega = cov_omega,
               # prop=prop, #proportion of intensity,
                #true_values = c(true1, true2),
                Validation_C = Validation_C,
                Validation_Y = Validation_Y,
                validation_indices_for_mismatch = validation_indices_for_mismatch,
                index_for_splitting = index_for_splitting,
                validation_mismatchC = validation_mismatchC,
               precipitation = precipitation,
               mean_temperature = mean_temperature,
               experience = experience,
               no_previous_obs = no_previous_obs,
               validation_correct_match = validation_correct_match,
               validation_correctmatchC = validation_correctmatchC,
               validation_correctmatchY = validation_correctmatchY)

distance_to_road <- list(training_data_with_distance, 
                         validation_data_with_distance)

#load("distance_to_road.RData")
save(distance_to_road, file="distance_to_road.RData")
save(dataout, file = "nimble_data.RData")
