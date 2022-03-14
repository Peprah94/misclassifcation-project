# Rscript that formats the data for analysis.
# Requires you to read in the data: viceroy.csv, queen.csv and monarch.csv

# Packages needed to run the script
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
viceroy <- read_csv("viceroy.csv")
monarch <- read_csv("monarch.csv")
queen <- read_csv("queen_new.csv")


#############
# LIMENITIS
#############
# Remove NAs and select data from the insect taxon
viceroy <- viceroy%>%
  filter(!is.na(taxon_genus_name))%>%
  filter(!is.na(taxon_species_name))%>%#remove NAs
  filter(iconic_taxon_name=="Insecta")#remove NAs

#setting the other genus as others
other <- filter(viceroy, !taxon_genus_name %in% c("Limenitis" ,
                                                  "Danaus"))$taxon_genus_name
viceroy[viceroy$taxon_genus_name %in%other[!is.na(other)],]$taxon_species_name <- "other"

#selecting the data that belongs to other limenitis
other_limenitis <- viceroy %>%
  filter(taxon_species_name != "Limenitis archippus")%>%
  filter(taxon_genus_name != "Danaus")%>%
  filter(taxon_species_name != "other")%>%
  dplyr::select(taxon_species_name)

viceroy[viceroy$taxon_species_name%in%(other_limenitis$taxon_species_name)[!is.na(other_limenitis$taxon_species_name)],]$taxon_species_name <- "other limenitis"

#Selecting the data that belongs to other danaus
other_danaus <- viceroy%>%
  filter(taxon_species_name == "Danaus eresimus")

viceroy[viceroy$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"

# Seperating the data into year, month and day
viceroy <- viceroy%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

# Assign a new column with the true name
viceroy$truth <- rep("Limenitis archippus", length(viceroy$taxon_species_name))


#############
# MONARCH
#############
# Remove NAs and select data from the insect taxon
monarch <- monarch%>%
  filter(!is.na(taxon_genus_name))%>%
  filter(!is.na(taxon_species_name))%>%#remove NAs
  filter(iconic_taxon_name=="Insecta")#remove NAs

#setting the other genus as others
other <- filter(monarch, !taxon_genus_name %in% c("Limenitis" ,
                                                  "Danaus"))$taxon_genus_name
monarch[monarch$taxon_genus_name %in%other[!is.na(other)],]$taxon_species_name <- "other"

#selecting the data that belongs to other limenitis
other_limenitis <- monarch %>%
  filter(taxon_species_name != "Limenitis archippus")%>%
  filter(taxon_genus_name != "Danaus")%>%
  filter(taxon_species_name != "other")%>%
  dplyr::select(taxon_species_name)

monarch[monarch$taxon_species_name%in%(other_limenitis$taxon_species_name)[!is.na(other_limenitis$taxon_species_name)],]$taxon_species_name <- "other limenitis"

#Selecting the data that belongs to other danaus
other_danaus <- monarch%>%
  filter(!taxon_species_name %in% c("other limenitis", "other"))%>%
  filter(taxon_genus_name != "Limenitis")%>%
  filter(!taxon_species_name %in% c("Danaus gilippus",  "Danaus plexippus"))

monarch[monarch$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"

# Seperating the data into year, month and day
monarch <- monarch%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

# Assign a new column with the true name
monarch$truth <- rep("Danaus plexippus", length(monarch$taxon_species_name))

#############
#QUEEN
#############
# Remove NAs and select data from the insect taxon
queen <- queen%>%
  filter(!is.na(taxon_genus_name))%>%
  filter(!is.na(taxon_species_name))%>%#remove NAs
  filter(iconic_taxon_name=="Insecta")

#setting the other genus as others
other <- filter(queen, !taxon_genus_name %in% c("Limenitis" ,
                                                "Danaus"))$taxon_genus_name
queen[queen$taxon_genus_name %in%other[!is.na(other)],]$taxon_species_name <- "other"

#selecting the data that belongs to other limenitis
other_limenitis <- queen %>%
  filter(taxon_species_name != "Limenitis archippus")%>%
  filter(taxon_genus_name != "Danaus")%>%
  filter(taxon_species_name != "other")%>%
  dplyr::select(taxon_species_name)

queen[queen$taxon_species_name%in%(other_limenitis$taxon_species_name)[!is.na(other_limenitis$taxon_species_name)],]$taxon_species_name <- "other limenitis"

#Selecting the data that belongs to other danaus
other_danaus <- queen%>%
  filter(!taxon_species_name %in% c("other limenitis", "other"))%>%
  filter(taxon_genus_name != "Limenitis")%>%
  filter(!taxon_species_name %in% c("Danaus gilippus",  "Danaus plexippus"))

queen[queen$taxon_species_name%in%(other_danaus$taxon_species_name)[!is.na(other_danaus$taxon_species_name)],]$taxon_species_name <- "other danaus"

# Seperating the data into year, month and day
queen <- queen%>%
  dplyr::arrange(id)%>%
  dplyr::mutate(year = lubridate::year(observed_on), 
                month = lubridate::month(observed_on), 
                day = lubridate::day(observed_on))

# Assign a new column with the true name
queen$truth <- rep("Danaus gilippus", length(queen$taxon_species_name))

#####################
#Putting all together
#######################
butterfly <- rbind(viceroy, monarch,queen)%>%
  filter(!is.na(longitude) |!is.na(longitude))

#function to subset data for each year
data_for_year <- function(year_id, butterfly){
  year_data <-butterfly%>%
    dplyr::filter(year==year_id)
  return(year_data=year_data)
}

#The years under consideration
unique_year <- c("2019", "2020")

# Generating the data
sim <- pblapply(unique_year, function(x){
  ret <- data_for_year(x,butterfly)
}, cl=1)

#save the data for use
save(sim, file="yearly_data.RData")
