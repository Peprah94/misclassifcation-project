source("nimble.R")

main <- run_simulations(dataout, "main", FALSE)

save(main, file="estimated_data_main1.RData")