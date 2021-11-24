set.seed(2021)
library(readr)
library(ggplot2)
library(dplyr)
library(pbapply)
library(spdep)
library(unix)
library("rnaturalearth")
library("rnaturalearthdata")
theme_set(theme_bw())
unix::rlimit_as(100*10^9)

titmouse<- read_csv("~/Documents/GitHub/Thesis/INaturalist/observations-150335.csv")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
titmouse<- read_csv("observations-150335.csv")

#Changing the names of some scientific names
titmouse$scientific_name[titmouse$scientific_name == "Baeolophus ridgwayi ridgwayi"] <- "Baeolophus ridgwayi"

#Unique species and scientific name
unique(titmouse$species_guess)
unique(titmouse$scientific_name)


#Select only those with the specified scientific name
titmouse <- titmouse%>%
  filter(scientific_name %in% c("Baeolophus bicolor","Baeolophus inornatus","Baeolophus atricristatus",
                                "Baeolophus wollweberi","Baeolophus ridgwayi"))

titmouse[grepl("uft",titmouse$species_guess),]$species_guess <- "tufted titmouse"
titmouse[grepl("bicolor",titmouse$species_guess),]$species_guess <- "tufted titmouse"
titmouse[grepl("ak ",titmouse$species_guess),]$species_guess <- "oak titmouse"
titmouse[grepl("AK ",titmouse$species_guess),]$species_guess <- "oak titmouse"
titmouse[grepl("inornatus",titmouse$species_guess),]$species_guess <- "oak titmouse"
titmouse[grepl("unip",titmouse$species_guess),]$species_guess <- "juniper titmouse"
titmouse[grepl("ridle",titmouse$species_guess),]$species_guess <- "bridled titmouse"
titmouse[grepl("wollweberi",titmouse$species_guess),]$species_guess <- "bridled titmouse"
titmouse[grepl("crested",titmouse$species_guess),]$species_guess <- "black-crested titmouse"
titmouse[grepl("atricrista",titmouse$species_guess),]$species_guess <- "black-crested titmouse"

other <- filter(titmouse, !species_guess %in% c("oak titmouse" ,
                                                "bridled titmouse",
                                                "juniper titmouse",
                                                "tufted titmouse" ,
                                                "black-crested titmouse" ))$species_guess
titmouse[titmouse$species_guess %in%other[!is.na(other)],]$species_guess <- "other"

library(tidyr)

#ind <- sample(seq(1:nrow(titmouse)),400, replace = FALSE)
# Data format
#Seperating according to years
titmouse <- titmouse %>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

year_data <-titmouse%>%
  dplyr::filter(year==2020)
#unique(titmouse$species_guess)
#titmouse%>%
#filter(grepl("ufted", species_guess))


#titmouse <- titmouse%>%
#filter(place_country_name =="United States")
par(mar = c(5.1, 5.1, 5.1, 5.1))
ggplot(data = world)+
  geom_sf()+
  geom_point(data =year_data[!is.na(year_data$species_guess),], aes(as.numeric(longitude), as.numeric(latitude), col=species_guess, shape=species_guess))+
  coord_sf(xlim = c(-125, -68), ylim = c(17.5, 46.9), expand = FALSE) +
  xlab("longitude")+
  ylab("latitude")#+
  theme(legend.position = c(0.8,0.2))#+
  facet_wrap(vars(scientific_name))
