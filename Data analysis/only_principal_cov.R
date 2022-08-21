source("nimble.R")

only_principal_cov <- run_simulations(dataout, "only_principal_cov", FALSE)

save(only_principal_cov , file="estimated_data_only_principal_cov1.RData")