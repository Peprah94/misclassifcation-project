source("nimble.R")

constant <- run_simulations(dataout, "constant", FALSE)

save(constant, file="estimated_data_constant1.RData")