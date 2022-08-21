load("simulated data1.RData")
source("estimation_simulation.R")

sim <- sim[1:200]

cl <- makeCluster(5)
clusterExport(cl, "run_simulations")
setDefaultCluster(cl)
rep_estimates <- pblapply(sim,function(x){
  ret <- run_simulations(x, type = "constant", model_selection = FALSE)
} , cl=cl)

save(rep_estimates, file="estimateddata_constant1.RData")
stopCluster(cl)