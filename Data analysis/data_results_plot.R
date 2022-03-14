# Script for plotting and summarizing results from the fitting_data_model.R
#plots percentage of correct classification 
#and the classification probabilities as a csv file

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

#world <- ne_countries(scale = "medium", returnclass = "sf")
#class(world)
load("needs_20.RData")
load("pred_needs_19.RData")
load("research_20.RData")
load("omega.RData")

#extracting data 
percentage_correct_needsid_19 <- pred_true_species_2019_needs_19$percentage_correct_after_pred
validation_percentage_correct_needsid_19 <- pred_true_species_2019_needs_19$percentage_correct_before_pred
percentage_correct_researchgrade_20 <- pred_true_species_2020_research_20$percentage_correct_after_pred
validation_percentage_correct_researchgrade_20 <- pred_true_species_2020_research_20$percentage_correct_before_pred
percentage_correct_needid_20 <- pred_true_species_2020_needs_20$percentage_correct_after_pred
validation_percentage_correct_needid_20 <- pred_true_species_2020_needs_20$percentage_correct_before_pred

#summarizing in a dataframe
all_percentages <- data.frame(data = rbind(percentage_correct_needsid_19,percentage_correct_researchgrade_20, percentage_correct_needid_20,
                                           validation_percentage_correct_needsid_19,validation_percentage_correct_researchgrade_20,validation_percentage_correct_needid_20),
                              year = as.factor(c("need id 2019", "research 2020", "need id 2020","need id 2019", "research 2020", "need id 2020")),
                              type = c("% after prediction","% after prediction" ,"% after prediction","% before prediction", "% before prediction","% before prediction"))%>%
  dplyr::mutate(data= data*100)

#plotting the results
ggplot(data=all_percentages, aes(x=year, y= data, fill= type))+
  geom_bar(stat = "identity", position=position_dodge2(reverse = TRUE))+
  geom_text(aes(label=round(data,2)), vjust=2.5, color="white",
            position = position_dodge2(1, reverse = TRUE), size=3)+
  scale_fill_brewer(palette="Paired", name="")+
  ylab("Percentage of correct classification")+
  xlab("Data type")

# Classification probabilites
row.names(omega) <- c("D. gilippus", "D. plexippus", "L. archipus")
colnames(omega) <- c("D. gilippus", "D. plexippus", "L. archipus", "other", "other Limenitis")
write.csv(omega, file="classification_probability.csv")

data_omega <- read.csv("classification_probability.csv")
colnames(data_omega) <- c("True state","D. gilippus", "D. plexippus", "L. archipus", "other", "other Limenitis")
write.csv(data_omega, file="classification_probability_data.csv")
