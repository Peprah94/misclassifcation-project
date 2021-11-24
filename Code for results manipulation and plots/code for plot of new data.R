library(readr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gganimate)
library(parallel)
library(doParallel)
library(pbapply)
library("rnaturalearth")
library("rnaturalearthdata")
theme_set(theme_bw())
unix::rlimit_as(100*10^9)


monthly_data <- function(i, plot = TRUE){
  library(readr)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
  library(gganimate)
  library(parallel)
  library(pbapply)
  library("rnaturalearth")
  library("rnaturalearthdata")
  theme_set(theme_bw())
  unix::rlimit_as(100*10^9)
  
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#read in data
viceroy <- read_csv("~/Documents/GitHub/Thesis/INaturalist/viceroy.csv")
monarch <- read_csv("~/Documents/GitHub/Thesis/INaturalist/monarch.csv")
queen <- read_csv("~/Documents/GitHub/Thesis/INaturalist/queen_new.csv")
#sa_monarch <- read_csv("~/Documents/GitHub/Thesis/INaturalist/South American Monarch.csv")


#############
# LIMENITIS
#############
#Formatting the data
viceroy <- viceroy%>%
  filter(!is.na(taxon_genus_name))%>%
  filter(!is.na(taxon_species_name))%>%#remove NAs
  filter(iconic_taxon_name=="Insecta")#remove NAs

#setting the other genus as others
other <- filter(viceroy, !taxon_genus_name %in% c("Limenitis" ,
                                                "Danaus"))$taxon_genus_name
viceroy[viceroy$taxon_genus_name %in%other[!is.na(other)],]$taxon_species_name <- "other"


other_limenitis <- viceroy %>%
  filter(taxon_species_name != "Limenitis archippus")%>%
  filter(taxon_genus_name != "Danaus")%>%
  filter(taxon_species_name != "other")%>%
  dplyr::select(taxon_species_name)

viceroy[viceroy$taxon_species_name%in%(other_limenitis$taxon_species_name)[!is.na(other_limenitis$taxon_species_name)],]$taxon_species_name <- "other limenitis"

other_danaus <- viceroy%>%
  filter(taxon_species_name == "Danaus eresimus")

viceroy[viceroy$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"

#Counting the number in each taxa
#tt <-viceroy%>%
#  count(taxon_species_name)

#Plotting of viceroy

viceroy <- viceroy%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

#year_data_viceroy <-viceroy%>%
 # dplyr::filter(year==2021)

#data_viceroy <- year_data_viceroy%>%
#  filter(longitude>= -120 &longitude <= -100)%>%
#  filter(latitude >= 25 & latitude <= 40)
viceroy$truth <- rep("Limentis archippus", length(viceroy$taxon_species_name))

#tt <-year_data_viceroy%>%
#  count(taxon_species_name)

#par(mar = c(5.1, 5.1, 5.1, 5.1))
#gg1 <- ggplot(data = world)+
 # geom_sf()+
 # geom_point(data =year_data_viceroy, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
 # coord_sf(xlim = c(-135, -50), ylim = c(17.5, 58), expand = FALSE) +
 # geom_rect(aes(xmin=-120, xmax=-100, ymin= 25, ymax=40), color="blue", alpha=0)+
 # xlab("longitude")+
 # ylab("latitud#e")

#############
# MONARCH
#############
monarch <- monarch%>%
  filter(!is.na(taxon_genus_name))%>%
  filter(!is.na(taxon_species_name))%>%#remove NAs
  filter(iconic_taxon_name=="Insecta")#remove NAs

#setting the other genus as others
other <- filter(monarch, !taxon_genus_name %in% c("Limenitis" ,
                                                  "Danaus"))$taxon_genus_name
monarch[monarch$taxon_genus_name %in%other[!is.na(other)],]$taxon_species_name <- "other"

other_limenitis <- monarch %>%
  filter(taxon_species_name != "Limenitis archippus")%>%
  filter(taxon_genus_name != "Danaus")%>%
  filter(taxon_species_name != "other")%>%
  dplyr::select(taxon_species_name)

monarch[monarch$taxon_species_name%in%(other_limenitis$taxon_species_name)[!is.na(other_limenitis$taxon_species_name)],]$taxon_species_name <- "other limenitis"

other_danaus <- monarch%>%
  filter(!taxon_species_name %in% c("other limenitis", "other"))%>%
  filter(taxon_genus_name != "Limenitis")%>%
  filter(!taxon_species_name %in% c("Danaus gilippus",  "Danaus plexippus"))

monarch[monarch$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"

monarch <- monarch%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

#year_data_monarch <-monarch%>%
#  dplyr::filter(year==2021)

#data_monarch <- year_data_monarch%>%
 # filter(longitude>= -120 &longitude <= -100)%>%
 # filter(latitude >= 25 & latitude <= 40)
monarch$truth <- rep("Danaus plexippus", length(monarch$taxon_species_name))

#tt <-year_data_monarch%>%
 # count(taxon_species_name)

#par(mar = c(5.1, 5.1, 5.1, 5.1))
#gg2 <-ggplot(data = world)+
 # geom_sf()+
 # geom_point(data =year_data_monarch, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
 # coord_sf(xlim = c(-135, -50), ylim = c(17.5, 58), expand = FALSE) +
 # geom_rect(aes(xmin=-120, xmax=-100, ymin= 25, ymax=40), color="blue", alpha=0)+
 # xlab("longitude")+
 # ylab("latitude")



#############
#QUEEN
#############
queen <- read_csv("~/Documents/GitHub/Thesis/INaturalist/queen.csv")
queen <- queen%>%
  filter(!is.na(taxon_genus_name))%>%
  filter(!is.na(taxon_species_name))%>%#remove NAs
  filter(iconic_taxon_name=="Insecta")
#setting the other genus as others
other <- filter(queen, !taxon_genus_name %in% c("Limenitis" ,
                                                     "Danaus"))$taxon_genus_name
queen[queen$taxon_genus_name %in%other[!is.na(other)],]$taxon_species_name <- "other"

other_limenitis <- queen %>%
  filter(taxon_species_name != "Limenitis archippus")%>%
  filter(taxon_genus_name != "Danaus")%>%
  filter(taxon_species_name != "other")%>%
  dplyr::select(taxon_species_name)

queen[queen$taxon_species_name%in%(other_limenitis$taxon_species_name)[!is.na(other_limenitis$taxon_species_name)],]$taxon_species_name <- "other limenitis"

other_danaus <- queen%>%
  filter(!taxon_species_name %in% c("other limenitis", "other"))%>%
  filter(taxon_genus_name != "Limenitis")%>%
  filter(!taxon_species_name %in% c("Danaus gilippus",  "Danaus plexippus"))

queen[queen$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"

queen <- queen%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

#year_data_queen <-queen%>%
 # dplyr::filter(year==2021)


#tt <-year_data_queen%>%
#  count(taxon_species_name)

#par(mar = c(5.1, 5.1, 5.1, 5.1))
#gg4 <- ggplot(data = world)+
 # geom_sf()+
 # geom_point(data =year_data_queen, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
 # coord_sf(xlim = c(-135, -50), ylim = c(17.5, 58), expand = FALSE) +
 # geom_rect(aes(xmin=-120, xmax=-100, ymin= 25, ymax=40), color="blue", alpha=0)+
 # xlab("longitude")+
 # ylab("latitude")

#ggpubr::ggarrange(gg1,gg2,gg3,gg4, common.legend = TRUE)

#data_queen <- year_data_queen%>%
 # filter(longitude>= -120 &longitude <= -100)%>%
 # filter(latitude >= 25 & latitude <= 40)

queen$truth <- rep("Danaus gilippus", length(queen$taxon_species_name))

all_data_together <- rbind(viceroy, monarch, queen)%>%
  filter(year%in%c("2019", "2020"))%>%
filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)

all_data <- rbind(viceroy, monarch, queen)%>%
  filter(year==2019)%>%
  filter(!is.na(longitude)|!is.na(latitude))%>%
  filter(quality_grade=="research")%>%
  filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)%>%
  filter(month==i)
  #filter(month %in% c("8"))
#usa=getData('GADM', country='USA', level=2)
#max_temp=getData('worldclim', var='tmax', res=10) 
#coords <- data.frame(cbind(all_data$longitude,all_data$latitude ))
#pts <- SpatialPointsDataFrame(coords = coords, data = all_data[,10:11], proj4string = crs(usa))
#dat <- intersect(pts@coords, usa@data)

#tt <-all_data%>%
 # count(taxon_species_name)

all_data_needs_id <- rbind(viceroy, monarch, queen)%>%
 filter(year==2019)%>%
filter(!is.na(longitude)|!is.na(latitude))%>%
filter(quality_grade=="needs_id")%>%
 filter(longitude < 60)%>%
 filter(latitude < 50 & latitude > 30)%>%
filter(month==i)

all_data_research_grade_19 <- rbind(viceroy, monarch, queen)%>%
  filter(year==2019)%>%
  filter(!is.na(longitude)|!is.na(latitude))%>%
  filter(quality_grade=="research")%>%
  filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)%>%
  filter(month==i)

all_data_research_grade_20 <- rbind(viceroy, monarch, queen)%>%
  filter(year==2020)%>%
  filter(!is.na(longitude)|!is.na(latitude))%>%
  filter(quality_grade=="research")%>%
  filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)%>%
  filter(month==i)

library(sp)
library(raster)
#init_plot <- SpatialPoints(coords = cbind( all_data$longitude,all_data$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#plot(init_plot)
#States <-raster::getData("GADM",country="United States",level=1)
#Florida <- States[States$NAME_1=="Florida",]
#points_in_norway <- raster::intersect(init_plot, Florida)
#coord_florida <- as.data.frame(points_in_norway@coords)

#all_data_florida <- all_data%>%
 # filter(longitude %in% coord_florida$coords.x1 & latitude %in% coord_florida$coords.x2)


#extract <- raster::extract(points_in_norway, shp, df=TRUE, cellnumbers=TRUE)
#all_data_us <- confusion_us <- confusion_florida <- a <- list()

#Finding out the cross tabulation of values
confusion_us <- table(all_data$truth, all_data$taxon_species_name)
#confusion_florida <- table(all_data_florida$truth, all_data_florida$taxon_species_name)

#Plot of data
#par(mar = c(10, 10, 10, 10))
if(plot==TRUE){
ggplot(data = world)+
geom_sf()+
 geom_point(data =all_data_together, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
  coord_sf(xlim = c(-135, -50), ylim = c(10.5, 58), expand = FALSE) +
  #geom_rect(aes(xmin=-100, xmax=-75, ymin= 23, ymax=38), color="blue", alpha=0)+
  xlab("longitude")+
  ylab("latitude")+
  facet_wrap(~truth, ncol=2)+
  theme(legend.position = c(0.8,0.2))+
labs(shape="Reported species")+
  labs(color="Reported species")
  #gganimate::transition_states(month, transition_length = 2, state_length = 1)+
  #ggtitle('Month {closest_state}')
}
  #labs(title = "Month: {frame_time}")+
  #gganimate::shadow_mark(alpha=0.3,size=0.5)

all_data_us <- all_data


#monthly_data(viceroy,monarch,queen,1)

#animate(a, nframes = 400, fps=5)

#b <- ggplot(data = world)+
#  geom_sf()+
#  geom_point(data =all_data_florida, aes(as.numeric(longitude), as.numeric(latitude), col=taxon_species_name, shape=taxon_species_name))+
 # coord_sf(xlim = c(-84, -78), ylim = c(24, 32), expand = FALSE) +
  #geom_rect(aes(xmin=-100, xmax=-75, ymin= 23, ymax=38), color="blue", alpha=0)+
  #xlab("longitude")+
  #ylab("latitude")+
  #facet_wrap(~truth, ncol=2)+
  #gganimate::transition_states(month, transition_length = 2, state_length = 1)+
  #ggtitle('Month {closest_state}')
  #gganimate::transition_time(month)+
  #labs(title = "Month: {frame_time}")+
  #gganimate::shadow_mark(alpha=0.3,size=0.5)
#Selecting 100 sites
#ind <- sample(1:length(all_data$id),300, replace=FALSE)
#all_data <- all_data[ind,]
library(raster)
if(dim(all_data)[1] !=0){
max_temp=getData('worldclim', var='tmax', res=10) 
altitude = getData('worldclim', var='alt', res=10) 
coords <- data.frame(cbind(all_data$longitude,all_data$latitude ))
pts <- SpatialPointsDataFrame(coords = coords, data = all_data[,10:11], proj4string = crs(max_temp))
pts1 <- SpatialPointsDataFrame(coords = coords, data = all_data[,10:11], proj4string = crs(altitude))
temp_data = raster::extract(max_temp, pts, df=T)[,(i+1)]
altitude_data = raster::extract(altitude, pts1, df=T)[,2]
}else{
  temp_data=0
  altitude_data=0
}
#plot(pts_data)

# Validation data
if(nrow(all_data_needs_id)!=0){
coords_validation <- data.frame(cbind(all_data_needs_id$longitude,all_data_needs_id$latitude ))
pts <- SpatialPointsDataFrame(coords = coords_validation, data = all_data_needs_id[,10:11], proj4string = crs(max_temp))
pts1 <- SpatialPointsDataFrame(coords = coords_validation, data = all_data_needs_id[,10:11], proj4string = crs(altitude))
temp_data_validation = raster::extract(max_temp, pts, df=T)[,(i+1)]
altitude_data_validation = raster::extract(altitude, pts1, df=T)[,2]
}else{
  temp_data_validation =0
  altitude_data_validation =0
  coords_validation=0
}
#Research grade for 2019
if(nrow(all_data_research_grade_19)!=0){
coords_19 <- data.frame(cbind(all_data_research_grade_19$longitude,all_data_research_grade_19$latitude ))
pts <- SpatialPointsDataFrame(coords = coords_19, data = all_data_research_grade_19[,10:11], proj4string = crs(max_temp))
pts1 <- SpatialPointsDataFrame(coords = coords_19, data = all_data_research_grade_19[,10:11], proj4string = crs(altitude))
temp_data_19 = raster::extract(max_temp, pts, df=T)[,(i+1)]
altitude_data_19 = raster::extract(altitude, pts1, df=T)[,2]
}else{
  temp_data_19 =0
  altitude_data_19 =0
  coords_19 =0
}
#Research grade for 2020
if(nrow(all_data_research_grade_20)!=0){
coords_20 <- data.frame(cbind(all_data_research_grade_20$longitude,all_data_research_grade_20$latitude ))
pts <- SpatialPointsDataFrame(coords = coords_20, data = all_data_research_grade_20[,10:11], proj4string = crs(max_temp))
temp_data_20 = raster::extract(max_temp, pts, df=T)[,(i+1)]
altitude_data_20 = raster::extract(altitude, pts, df=T)[,2]
}else{
  temp_data_20 =0
  altitude_data_20 =0
  coords_20 =0
}

### Scaling of temp data
mean_temp_19 <- mean(temp_data, na.rm=TRUE)
mean_alt_19 <- mean(altitude_data,na.rm=TRUE)
sd_temp_19 <- sd(temp_data, na.rm = TRUE)
sd_alt_19 <- sd(altitude_data, na.rm = TRUE)
temp = (temp_data-mean_temp_19)/sd_temp_19
altitude = (altitude_data-mean_alt_19)/sd_alt_19


####

returned_dataframe <- data.frame(C=as.numeric(as.factor(all_data$truth)), 
                                 Y=as.numeric(as.factor(all_data$taxon_species_name)),
                                 agree = all_data$num_identification_agreements,
                                 temp = temp,
                                 altitude = altitude)%>%
  filter(!is.na(temp) | !is.na(altitude))

#locs <- cbind(all_data_florida[[1]]$latitude,all_data_florida[[1]]$longitude)
data_for_nimble <- list(C=returned_dataframe$C,
                        Y = returned_dataframe$Y,
                        #locs=locs,
                        confusion_us = confusion_us,
                        #cov=as.numeric(as.factor(all_data_florida$quality_grade)),
                        agree= returned_dataframe$agree,
                        disagree= all_data$num_identification_disagreements,
                        temp = returned_dataframe$temp,
                        altitude = returned_dataframe$altitude,
                        validation_C = as.numeric(as.factor(all_data_needs_id$truth)),
                        validation_Y = as.numeric(as.factor(all_data_needs_id$taxon_species_name)),
                        validation_agree = all_data_needs_id$num_identification_agreements,
                        validation_disagree= all_data_needs_id$num_identification_disagreements,
                        temp_validation = (temp_data_validation-mean_temp_19)/sd_temp_19,
                        altitude_validation =(altitude_data_validation-mean_alt_19)/sd_alt_19,
                        coords_validation=coords_validation,
                        research_C19 = as.numeric(as.factor(all_data_research_grade_19$truth)),
                        research_Y19 = as.numeric(as.factor(all_data_research_grade_19$taxon_species_name)),
                        research19_agree = all_data_research_grade_19$num_identification_agreements,
                        research19_disagree= all_data_research_grade_19$num_identification_disagreements,
                        temp_19 = (temp_data_19-mean_temp_19)/sd_temp_19,
                        altitude_19 =(altitude_data_19-mean_alt_19)/sd_alt_19,
                        coords_19=coords_19 ,
                        research_C20 = as.numeric(as.factor(all_data_research_grade_20$truth)),
                        research_Y20 = as.numeric(as.factor(all_data_research_grade_20$taxon_species_name)),
                        research20_agree = all_data_research_grade_20$num_identification_agreements,
                        research20_disagree= all_data_research_grade_20$num_identification_disagreements,
                        temp_20 = (temp_data_20-mean_temp_19)/sd_temp_19,
                        altitude_20 = (altitude_data_20-mean_alt_19)/sd_alt_19,
                        coords_20=coords_20 )
return(data_for_nimble)
}



months <- as.list(seq(1,12,1))
cl <- makeCluster(2)
setDefaultCluster(cl)
data_for_nimble <- pblapply(months,monthly_data, cl=cl)
save(data_for_nimble, file="data_new.RData")
stopCluster(cl)


