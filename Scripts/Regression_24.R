rm(list=ls())
source("Scripts/Utility_Functions.R")

# Regression models on Provinces
# Create the regression models for each year and class. Set some parameters and utilities:

classe = c("0_10","11_14","15_21")
years = paste0("T_",11:24)
knots = 52
xcp = as.numeric((as.numeric(as.Date("2020-01-01")):as.numeric(as.Date("2020-12-31"))))
knots_vector = floor(xcp[seq(1,length(xcp), length.out = knots)])
to_average = 4
alpha.seq = seq(0, 1, by=0.05)  ##alpha transformation parameter
class = classe[3]
window = 1:158
years_names = paste0(20, substr(years,3,4))


##load past years best models
load("Output/Data/Regression_Province_clr_23.Rdata")
load("Output/Data/Regression_Province_best_alpha_23.Rdata")
load("Output/Data/Regression_Province_best_alpha_isometric_23.Rdata")


#######FUNCTION

models.predictions = function(alpha, transformation = "Tsagris") {
  idx = paste0(alpha)

  i = length(years)
  
  year = years[i]
  argvals = clist_province[[idx]][[years[i-1]]]$argvals
  argvals_restricted = clist_province[[idx]][[year]]$argvals
  
  mean_previous_provinces =  clist_province[[idx]][[years[i-1]]]$data
  for (y in years[(i - to_average):(i-2)]) 
    mean_previous_provinces = clist_province[[idx]][[y]]$data + mean_previous_provinces
  mean_previous_provinces = mean_previous_provinces/to_average
  mean_previous_provinces = fdata(mean_previous_provinces, argvals)
  
  model = regr_models_provinces_best[[years[i-1]]]$model
  if (transformation == "Isometric") model = regr_models_provinces_best_isom[[years[i-1]]]$model
  if (transformation == "CLR") model = regr_models_provinces_clr[[years[i-1]]]$model
  
  temp = predict(model, mean_previous_provinces)
  preds = fdata(temp$data[,window], argvals_restricted)
  
  errors = clist_province[[idx]][[year]] - preds
  
  meape = apply(abs(errors$data)/abs(clist_province[[idx]][[year]]$data), 1, median)
  
  preds_density = NULL
  if (transformation == "Tsagris") preds_density = alpha.folding(preds, alpha)
  if (transformation == "Isometric") preds_density = alpha.isometric.inv(preds, alpha)
  if (transformation == "CLR") preds_density = inv_clr(preds)
  
  return(list(predictions = preds, errors = errors, MedAPE = meape))
}


####### TSAGRIS

load("Output/Data/Smoothing_Provinces_8_alpha_24.Rdata")  #smoothed curves

tsagris_performance = models.predictions(alpha.best, "Tsagris")

####### CLR
clr_performance = models.predictions(0, "CLR")

####### ISOMETRIC
load("Output/Data/Smoothing_Provinces_8_alpha_isometric_24.Rdata")

isometric_performance = models.predictions(alpha.best_isom, "Isometric")

save(clr_performance, tsagris_performance, isometric_performance, file = Rdata_path("Predictive_performance_24"))

######## PLOTS

rm(list=ls())
load("Output/Data/Predictive_performance_24.Rdata")
source("Scripts/Utility_Functions.R")


####### Density comparison plots for KLdiv and MedAPE

constant_density = fdata(rep(1/diff(clr_performance$predictions$rangeval), 
                             length(clr_performance$predictions$argvals)),
                         clr_performance$predictions$argvals)

KL_div = function(errors, baseline_density, transformation, alpha = 0) {
  errors_density = NULL
  
  if (transformation == "CLR") errors_density = inv_clr(errors)
  if (transformation == "Tsagris") errors_density = alpha.folding(errors, alpha)
  if (transformation == "Isometric") errors_density = alpha.isometric.inv(errors, alpha)
  
  return( int.simpson(errors_density * log(errors_density/baseline_density)))
  
}

years_plot = "T_24"
years_plot_name = "2024"
comparison = data.frame(Year = rep(years_plot_name, each = 107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), each=107),
                        KLdiv = NA, 
                        MedAPE = NA)

alpha.best = 0.85
alpha.best_isom = 1

comparison[comparison$Transformation == "CLR" & comparison$Year == years_plot_name,3] = 
  KL_div(clr_performance$errors, constant_density, "CLR")

comparison[comparison$Transformation == "$A_{0.85}$" & comparison$Year == years_plot_name,3] =    
  KL_div(tsagris_performance$errors, constant_density, "Tsagris", alpha.best)

comparison[comparison$Transformation == "$A_{1-IT}$" & comparison$Year == years_plot_name,3] = 
  KL_div(isometric_performance$errors, constant_density, "Isometric", alpha.best_isom)


g1 =ggplot(comparison, aes(x = as.factor(Year), y = KLdiv, fill = Transformation)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Align the boxplots
  labs(x = "", y = "KL divergence") +
  scale_fill_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                    labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$")))  +
  theme_minimal() + ylim(c(0,0.06)) + 
  theme(legend.position = "none") 


comparison[comparison$Transformation == "CLR",4] = clr_performance$MedAPE
comparison[comparison$Transformation == "$A_{0.85}$",4] = tsagris_performance$MedAPE
comparison[comparison$Transformation == "$A_{1-IT}$",4] = isometric_performance$MedAPE

q1 =  ggplot(comparison, aes(x = as.factor(Year), y = MedAPE, fill = Transformation)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Align the boxplots
  labs(x = "", y = "MedAPE") +
  scale_fill_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                    labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  theme_minimal() + ylim(c(0,5)) + 
  theme(legend.position = "none") 

dev.new(width=15, height=7)
(g1 + q1 ) + plot_layout( guides = 'collect') & theme(legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_blank())



  

