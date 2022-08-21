
README.txt

A: Downloaded data from INaturalist

The data provided was downloaded on 28.10.2021 from INaturalist using the query below

a. Monarch.csv
Query q=monarch&quality_grade=any&identifications=any&place_id=any&verifiable=true
Columns id, observed_on_string, observed_on, time_observed_at, time_zone, user_id, user_login, created_at, updated_at, quality_grade, license, url, image_url, sound_url, tag_list, description, num_identification_agreements, num_identification_disagreements, captive_cultivated, oauth_application_id, place_guess, latitude, longitude, positional_accuracy, public_positional_accuracy, geoprivacy, taxon_geoprivacy, coordinates_obscured, positioning_method, positioning_device, species_guess, scientific_name, common_name, iconic_taxon_name, taxon_id, taxon_genus_name, taxon_species_name

b. viceroy.csv
Query q=viceroy&quality_grade=any&identifications=any&place_id=any&verifiable=true
Columns id, observed_on_string, observed_on, time_observed_at, time_zone, user_id, user_login, created_at, updated_at, quality_grade, license, url, image_url, sound_url, tag_list, description, num_identification_agreements, num_identification_disagreements, captive_cultivated, oauth_application_id, place_guess, latitude, longitude, positional_accuracy, public_positional_accuracy, geoprivacy, taxon_geoprivacy, coordinates_obscured, positioning_method, positioning_device, species_guess, scientific_name, common_name, iconic_taxon_name, taxon_id, taxon_genus_name, taxon_species_name

c. queen.csv
Query q=Danaus+gilippus&search_on=names&quality_grade=any&identifications=any&place_id=any&verifiable=true
Columns id, observed_on_string, observed_on, time_observed_at, time_zone, user_id, user_login, created_at, updated_at, quality_grade, license, url, image_url, sound_url, tag_list, description, num_identification_agreements, num_identification_disagreements, captive_cultivated, oauth_application_id, place_guess, latitude, longitude, positional_accuracy, public_positional_accuracy, geoprivacy, taxon_geoprivacy, coordinates_obscured, positioning_method, positioning_device, species_guess, scientific_name, common_name, iconic_taxon_name, taxon_id, taxon_genus_name, taxon_species_name

B: R scripts
a. Simulation Study
The code should be run in this order for the results presented in the paper:
- simulation.R is used to simulate the data that is needed for the model.
- estimate_simulation.R is used to fit the model and estimate the parameters in NIMBLE
- plot_simulations.R is used to summarize the results from the model fitting in tables and graphs

b. Application to Butterfly data
The code should be run in this order for the results presented in the paper:
- Data_formatting.R is used to format the data that is needed for the model. Needs the queen.csv, viceroy.csv and monarch.csv as inputs.
- nimble.R is used to fit the model with INLA and NIMBLE.
- plot for data.R is used to summarise the results from nimble.R and present results in Tables and graphs.



