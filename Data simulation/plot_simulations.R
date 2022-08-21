library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 10))




niter = 200
load("simulated data1.RData")
load("estimateddata_variable1.RData")
variable <- rep_estimates
load("estimateddata_constant1.RData")
constant <- rep_estimates
load("estimateddata_intercept1.RData")
intercept <- rep_estimates
load("estimateddata_main1.RData")
main <- rep_estimates
load("estimateddata_only_principal_cov1.RData")
only_principal <- rep_estimates
load("estimateddata_opcovfinter1.RData")
opcovfinter <- rep_estimates

#Put all datasets together
all_data <- c(variable, intercept, constant, main,only_principal, opcovfinter)

#Functions to select parameters of interest from the simulation study
#For constant returned values
subset_parameters <- function(data, subset_index){
 ret <-  lapply(data, function(x){
    d <- x[[subset_index]]
  })
  return(ret)
}

#For constant returned values
subset_parameters_annex <- function(data, subset_index){
  ret <-  lapply(data, function(x){
    d <- c(x[[subset_index]])
  })
  return(ret)
}

#return accuracy
accuracy <- subset_parameters(all_data, 6)%>%
  unlist()

#return precision estimates
precision<- subset_parameters(all_data, 7)%>%
  unlist()

plot_data <- base::data.frame(accuracy=accuracy, precision=precision)
colnames(plot_data) <- c("accuracy", "precision")
plot_data1 <- plot_data%>%
  mutate(method = c(rep("CCM", niter), 
                    rep("CIM", niter),
                    rep("COM", niter),
                    rep("MCM", niter),
                    rep("CSCM", niter),
                    rep("CSICM", niter)),
         group = rep(c(rep("full", (niter/2)), 
                       rep("reduced", (niter/2))), 6))%>%
  reshape2::melt(id.vars=c("method", "group"))%>%
  filter(method != "MCM")%>%
  ggplot()+
  geom_boxplot(mapping = aes(x = as.factor(method), 
                             y = value, 
                             #col = as.factor(variable),
                             fill = as.factor(variable)))+
  ylab("Predictive value")+
  xlab("method")+
  theme(legend.title =element_blank(), legend.position = "top")+
  ylim(c(0.75,1))+
  #scale_y_cont+inuous(trans = "sqrt")+
  facet_wrap(~group)
ggsave("accuracy_and_precision.png", plot = plot_data1, dpi = 320, 
       height = 10, width = 10, units = "cm")

p <- subset_parameters(all_data, 8)%>%
  unlist()%>%
  data.frame()%>%
  mutate(method = c(rep("CCM", niter), rep("CIM", niter),rep("COM", niter),rep("MCM", niter),
                    rep("CSCM", niter),
                    rep("CSICM", niter)),
         group = rep(c(rep("full", (niter/2)), rep("reduced", (niter/2))), 6))%>%
  ggplot()+
  geom_boxplot(mapping = aes(x = as.factor(method), y = ., col = as.factor(method)))+
  ylab("Probability of selecting \n the classification covariate")+
  xlab("method")+
  theme(legend.position="none")+
  facet_wrap(~group)



coverage <- subset_parameters(all_data, 5)%>%
  do.call("rbind", .)%>%
  data.frame()%>%
  mutate(method = c(rep("CCM", niter), rep("CIM", niter),rep("COM", niter),rep("MCM", niter),
                    rep("CSCM", niter),
                    rep("CSICM", niter)),
         group = rep(c(rep("full", (niter/2)), rep("reduced", (niter/2))), 6))%>%
  group_by(group, method)%>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))%>%
  reshape2::melt(id.vars = c("group", "method"))%>%
  ungroup()%>%
  dplyr::filter(!variable %in% c("N_over_rhat", "all_rhat"))%>%
  mutate(variable = c(rep("beta_01", 12), rep("beta_11", 12), 
                      rep("omega0_11", 12), rep("omega0_21", 12),
                      rep("omega0_12", 12), rep("omega0_22", 12), 
                      rep("omega1_11", 12), rep("omega1_21", 12),
                      rep("omega1_12", 12), rep("omega1_22", 12)))%>%
  ggplot()+
  geom_boxplot(mapping = aes(x = as.factor(variable), y = value, col = as.factor(method)))+
  geom_hline(yintercept = 0.95, col = "black")+
  ylab("Coverage")+
  scale_x_discrete("Parameters", 
                     labels =  c("beta_01" = expression("\u03b2"[0]),
                                 "beta_11" = expression("\u03b2"[1]),
                                 "omega0_11" = expression("\u03c9"[0011]),
                                 "omega0_21" = expression("\u03c9"[0021]),
                                 "omega0_12" = expression("\u03c9"[0012]),
                                 "omega0_22" = expression("\u03c9"[0022]),
                                 "omega1_11" = expression("\u03c9"[111]),
                                 "omega1_21" = expression("\u03c9"[121]),
                                 "omega1_12" = expression("\u03c9"[112]),
                                 "omega1_22" = expression("\u03c9"[122])
                   
))+
 # xlab("Parameters")+
  facet_wrap(~group)
  

beta0 <- unlist(subset_parameters(all_data, 1))+ 1
beta1 <- unlist(subset_parameters(all_data, 2)) - 8
omega0 <- subset_parameters_annex(all_data, 3)%>%
  do.call("rbind", .)%>%
  data.frame()%>%
  mutate(X1 = X1 - sim[[1]]$true_values[1],
         X2 = X2 - sim[[1]]$true_values[2],
         X3 = X3 - sim[[1]]$true_values[3],
         X4 = X4 - sim[[1]]$true_values[4])

omega1 <- subset_parameters_annex(all_data, 4)%>%
  do.call("rbind", .)%>%
  data.frame()%>%
  mutate(X1 = X1 - sim[[1]]$true_values[5],
         X2 = X2 - sim[[1]]$true_values[6],
         X3 = X3 - sim[[1]]$true_values[7],
         X4 = X4 - sim[[1]]$true_values[8])

method = c(rep("CCM", niter), rep("CIM", niter),rep("COM", niter),rep("MCM", niter),
           rep("CSCM", niter),
           rep("CSICM", niter))
group = rep(c(rep("full", (niter/2)), rep("reduced", (niter/2))), 6)
bias_data <- cbind(beta0, beta1, omega0, omega1, method, group)
colnames(bias_data) <- c("beta_01", "beta_11", 
                         "omega0_11", "omega0_21",
                         "omega0_12", "omega0_22", 
                         "omega1_11", "omega1_21",
                         "omega1_12", "omega1_22", 
                         "method", "group")

#Setting some subsets of the data as NAs
bias_data1 <- bias_data%>%
  mutate(omega0_11 = ifelse(method %in% c("COM","MCM") & group %in% c("reduced","full"), NA, omega0_21),
         omega0_21 = ifelse(method %in% c("COM","MCM") & group %in% c("reduced","full"), NA, omega0_21),
         omega0_12 = ifelse(method %in% c("COM","MCM") & group %in% c("reduced","full"), NA, omega0_12),
         omega0_22 = ifelse(method %in% c("COM","MCM") & group %in% c("reduced","full"), NA, omega0_22),
         omega1_11 = ifelse(method %in% c("COM","MCM", "CIM") & group %in% c("reduced","full"), NA, omega1_11),
         omega1_12 = ifelse(method %in% c("COM","MCM", "CIM") & group %in% c("reduced","full"), NA, omega1_12),
         omega1_21 = ifelse(method %in% c("COM","MCM", "CIM") & group %in% c("reduced","full"), NA, omega1_21),
         omega1_22 = ifelse(method %in% c("COM","MCM", "CIM") & group %in% c("reduced","full"), NA, omega1_22),)



bias_plot <- bias_data1%>%
  reshape2::melt(id.vars = c("group", "method"))%>%
  ggplot()+
  geom_boxplot(mapping = aes(x = as.factor(variable), y = value, fill = as.factor(method)))+
  geom_hline(yintercept = 0, col = "black", linetype="dashed")+
  ylab("Bias")+
  scale_fill_discrete(name="method")+
  scale_x_discrete("Parameters", 
                   labels =  c("beta_01" = expression("\u03b2"[0]),
                               "beta_11" = expression("\u03b2"[1]),
                               "omega0_11" = expression("\u03c9"[0011]),
                               "omega0_21" = expression("\u03c9"[0021]),
                               "omega0_12" = expression("\u03c9"[0012]),
                               "omega0_22" = expression("\u03c9"[0022]),
                               "omega1_11" = expression("\u03c9"[111]),
                               "omega1_21" = expression("\u03c9"[121]),
                               "omega1_12" = expression("\u03c9"[112]),
                               "omega1_22" = expression("\u03c9"[122])
                               
                   ))+
  #scale_y_continuous(trans='log10')+
  # xlab("Parameters")+
  facet_wrap(~group)

main_fig_sim <- ggpubr::ggarrange(bias_plot, ncol=1, nrow = 2,
                  common.legend = TRUE)

ggsave("bias_plot_and_p.png", plot =main_fig_sim, dpi = 320, 
       height = 15, width = 12, units = "cm")

supp_fig_sim <- ggpubr::ggarrange(coverage, ncol=1, nrow = 2,
                  common.legend = TRUE)
ggsave("coverage.png", plot = supp_fig_sim, dpi = 320, 
       height = 10, width = 10, units = "cm")
