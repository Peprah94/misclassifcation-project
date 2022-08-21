source("nimble.R")

intercept <- run_simulations(dataout, "intercept", FALSE)

save(intercept, file="estimated_data_intercept1.RData")