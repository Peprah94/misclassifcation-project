source("nimble.R")

op_cov_finter <- run_simulations(dataout, "op_cov_finter", FALSE)

save(op_cov_finter, file="estimated_data_op_cov_finter1.RData")