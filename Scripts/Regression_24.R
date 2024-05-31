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


############################  comparison between CLR, Tsagris and Isometric

comparison = data.frame(Year = rep("2024", each=107), 
                        Transformation = rep(c("CLR", "Tsagris","Isometric"), each=107),
                        MedAPE = NA)

comparison[comparison$Transformation == "CLR",3] = clr_performance$MedAPE
comparison[comparison$Transformation == "Tsagris",3] = tsagris_performance$MedAPE
comparison[comparison$Transformation == "Isometric",3] = isometric_performance$MedAPE


ggplot(comparison, aes(x = as.factor(Year), y = MedAPE, fill = Transformation)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +  # Align the boxplots
  labs(x = "Year", y = "Value", title = 
         TeX(paste0("MedAPE comparison. Tsagris $\\alpha = 0.85$, Isometric $\\alpha = 1$"))) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("CLR" = "skyblue", "Tsagris" = "orange", "Isometric" = "red"))+
  ylim(c(0.5,6))


comparison = data.frame(
  Year = rep("2024", each=107), 
  Transformation = rep(c("CLR", "Tsagris","Isometric"), each=107),
  Kldiv = NA)


constant_density = fdata(rep(1/diff(clist_province[["0"]][["T_24"]]$rangeval), 
                             length(clist_province[["0"]][["T_24"]]$argvals)),
                         clist_province[["0"]][["T_24"]]$argvals)

KL_div = function(errors, baseline_density, transformation, alpha = 0) {
  errors_density = NULL
  
  if (transformation == "CLR") errors_density = inv_clr(errors)
  if (transformation == "Tsagris") errors_density = alpha.folding(errors, alpha)
  if (transformation == "Isometric") errors_density = alpha.isometric.inv(errors, alpha)
  
  return( int.simpson(errors_density * log(errors_density/baseline_density)))
  
}


comparison[comparison$Transformation == "CLR",3] = 
  KL_div(clr_performance$errors, constant_density, "CLR")

comparison[comparison$Transformation == "Tsagris",3] =    
  KL_div(tsagris_performance$errors, constant_density, "Tsagris", alpha.best)

comparison[comparison$Transformation == "Isometric",3] = 
  KL_div(isometric_performance$errors, constant_density, "Isometric", alpha.best_isom)



ggplot(comparison, aes(x = as.factor(Year), y = Kldiv, fill = Transformation)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +  # Align the boxplots
  #geom_vline(xintercept = seq(1.5, length(unique(comparison$Year)) - 0.5, by = 1), linetype = "dashed", color = "grey") +
  labs(x = "Year", y = "Value", title = 
         TeX(paste0("KL-divergence comparison. Tsagris $\\alpha = 0.85$, Isometric $\\alpha = 1$"))) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  # ylim(c(0,0.15)) +
  scale_fill_manual(values = c("CLR" = "skyblue", "Tsagris" = "orange", "Isometric" = "red"))



save(clr_performance, tsagris_performance, isometric_performance, file = Rdata_path("Predictive_performance_24"))

