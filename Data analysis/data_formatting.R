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
library(spatstat)
library(dplyr)
library(sp)
library(gstat)
library(ggplot2)
library(INLA)
library(rgeos)
theme_set(theme_bw())
unix::rlimit_as(100*10^9)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#read in data
viceroy <- read_csv("viceroy.csv")
monarch <- read_csv("monarch.csv")
queen <- read_csv("queen_new.csv")

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

viceroy <- viceroy%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))


viceroy$truth <- rep("Limenitis archippus", length(viceroy$taxon_species_name))


#############
# MONARCH
#############
monarch <- monarch%>%
  filter(!is.na(taxon_genus_name))%>%
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


monarch$truth <- rep("Danaus plexippus", length(monarch$taxon_species_name))



#############
#QUEEN
#############
queen <- read_csv("queen.csv")
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


queen$truth <- rep("Danaus gilippus", length(queen$taxon_species_name))

titmouse <- rbind(viceroy, monarch,queen)%>%
  filter(!is.na(longitude) |!is.na(longitude))%>%
  dplyr::group_by(user_id)%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(count = n(),
                previous_obs = row_number() - 1)%>%
  dplyr::ungroup()%>%
  dplyr::filter(year %in% c("2015", "2016","2017", "2018", "2019", "2020"))

gg1 <- ggplot(data = titmouse)+
  geom_bar(aes(x= previous_obs))+
  xlab("Number of previous reports")+
  ggtitle("a) Previous number of reports \n by citizen scientists")

gg2 <- ggplot(data = titmouse)+
  geom_bar(aes(x= count))+
  xlab("Total number of reports from each user")+
  ggtitle("b) Total number of reports \n by citizen scienti from 2017 to 2020")

cnt <- ggpubr::ggarrange(gg1, gg2, nrow=1, ncol=2)
ggsave("count_distribution.png", plot = cnt, dpi = 320, 
       height = 8, width = 10, units = "cm")
#Add the states to data
#state_of_interest = "Florida"
States <- raster::getData("GADM", country = "United States", level = 2)
#States <- States[ States$NAME_1 == state_of_interest,]
USborder <- rgeos::gUnaryUnion(States, id = States$ISO)
boundary_points <- fortify(USborder)
sp_df <- SpatialPointsDataFrame(cbind(as.numeric(titmouse$longitude), 
                                      as.numeric(titmouse$latitude)), 
                                titmouse,proj4string = CRS("+proj=longlat +datum=WGS84"))
df_titmouse <- raster::intersect(sp_df, States)
titmouse <- df_titmouse@data

alt<- raster::getData('worldclim', var='alt', res=10) 
bio <- raster::getData('worldclim', var='bio', res=10) 
sp_df <- SpatialPointsDataFrame(cbind(as.numeric(titmouse$longitude), 
                                      as.numeric(titmouse$latitude)), 
                                titmouse,proj4string = CRS("+proj=longlat +datum=WGS84"))
altitude <- raster::extract(alt, sp_df, df=T)
bio_data <- raster::extract(bio, sp_df, df = T)

titmouse <- cbind(titmouse, altitude, bio_data)%>%
  data.frame()


data_for_year <- function(year_id, titmouse){
  year_data <-titmouse%>%
    dplyr::filter(year==year_id)
  
  return(year_data=year_data)
}

unique_year <- c("2015", "2016", "2017", "2018", "2019")

training_data <- pblapply(unique_year, function(x){
  ret <- data_for_year(x,titmouse)
}, cl=1)%>%
  do.call("rbind", .)
#table(training_data$truth, training_data$taxon_species_name)
validation_data <- pblapply(list("2020"), function(x){
  ret <- data_for_year(x,titmouse)
}, cl=1)%>%
  do.call("rbind", .)


sim <- list(training_data, validation_data)

save(sim, file="yearly_data_new.RData")
