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

#
#Load results from data manipulation
full_data_result_19 <- read_csv("results_from_data.csv")
full_data_result_20 <- read_csv("results_from_data_2020.csv")

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
viceroy$truth <- rep("Limenitis archippus", length(viceroy$taxon_species_name))

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
  filter(latitude < 50 & latitude > 30)

all_data_needs_id <- rbind(viceroy, monarch, queen)%>%
  filter(year==2019)%>%
  filter(!is.na(longitude)|!is.na(latitude))%>%
  filter(quality_grade=="needs_id")%>%
  filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)

all_data_research_grade_19 <- rbind(viceroy, monarch, queen)%>%
  filter(year==2019)%>%
  filter(!is.na(longitude)|!is.na(latitude))%>%
  filter(quality_grade=="research")%>%
  filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)

all_data_research_grade_20 <- rbind(viceroy, monarch, queen)%>%
  filter(year==2020)%>%
  filter(!is.na(longitude)|!is.na(latitude))%>%
  filter(quality_grade=="research")%>%
  filter(longitude < 60)%>%
  filter(latitude < 50 & latitude > 30)


percentage_correct_needsid_19 <- sum(all_data_needs_id$truth==all_data_needs_id$taxon_species_name)/nrow(all_data_needs_id)
percentage_correct_researchgrade_20 <- sum(all_data_research_grade_20$truth == all_data_research_grade_20$taxon_species_name)/nrow(all_data_research_grade_20)


validation_percentage_correct_needsid_19 <-sum(full_data_result_19$truth==full_data_result_19$`max species`)/nrow(full_data_result_19) #proportion of max species that correspond to the truth



#monthly percentage of correct identifications
validation_prob <- full_data_result_19 %>%
  select(truth, `max species`, month) %>%
  dplyr:: group_by(month) %>%
  summarise(probability = mean(truth==`max species`, na.rm=TRUE) )


#percentage of correct identifications
validation_percentage_correct_needsid_20 <- sum(full_data_result_20$truth==full_data_result_20$max_species)/nrow(full_data)

#monthly percentage of correct identifications
research_2020 <-full_data_result_20 %>%
  select(truth, max_species, month) %>%
  dplyr:: group_by(month) %>%
  summarise(probability = mean(truth==max_species, na.rm=TRUE) )

all_percentages <- data.frame(data = rbind(percentage_correct_needsid_19,percentage_correct_researchgrade_20, validation_percentage_correct_needsid_19,validation_percentage_correct_needsid_20),
                              year = as.factor(c(2019, 2020, 2019, 2020)),
                              type = c("% before prediction", "% before prediction","% after prediction","% after prediction" ))%>%
  mutate(data= data*100)

ggplot(data=all_percentages, aes(x=year, y= data, fill= type))+
  geom_bar(stat = "identity", position=position_dodge2(reverse = TRUE))+
  geom_text(aes(label=round(data,2)), vjust=2.5, color="white",
            position = position_dodge2(1, reverse = TRUE), size=6)+
  scale_fill_brewer(palette="Paired", name="")+
  ylab("Percentage of correct classification")
