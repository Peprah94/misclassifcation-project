source("nimble.R")

variable <- run_simulations(dataout, "variable", FALSE)

save(variable, file="estimated_data_variable1.RData")