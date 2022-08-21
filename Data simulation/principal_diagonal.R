load("simulated data1.RData")
source("estimation_simulation.R")

sim <- sim[1:200]

cl <- makeCluster(5)
clusterExport(cl, "run_simulations")
setDefaultCluster(cl)
rep_estimates <- pblapply(sim,function(x){
  ret <- run_simulations(x, type = "only_principal_cov", model_selection = FALSE)
} , cl=cl)

save(rep_estimates, file="estimateddata_only_principal_cov1.RData")
stopCluster(cl)