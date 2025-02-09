rm(list=ls())
source("Scripts/Utility_Functions.R")

classe = c("0_10","11_14","15_21")
years = paste0("T_",11:23)
years_names = paste0(20,11:23)
knots = 52
xcp = as.numeric((as.numeric(as.Date("2020-01-01")):as.numeric(as.Date("2020-12-31"))))
knots_vector = floor(xcp[seq(1,length(xcp), length.out = knots)])
to_average = 4
alpha.seq = seq(0, 1, by=0.05)  ##alpha transformation parameter
class = classe[3]

ncores <- max(parallel::detectCores() -1,1)
if (ncores==1) {
  foreach::registerDoSEQ()
} else{
  cl <- suppressWarnings(parallel::makePSOCKcluster(ncores ))
  doParallel::registerDoParallel(cl)
}
ops.fda.usc(reset=TRUE, verbose=TRUE)

### FUNCTION TO FIT THE MODELS

model.fitting = function(alpha, transformation = "Tsagris") {
  idx = paste(alpha,"", sep="")
  
  regr_models_provinces = r.sq = meape =  list()
  
  for(i in (to_average+1):length(years))
  {
    year = years[i]
    argvals = clist_province[[idx]][[year]]$argvals
    mean_previous_provinces =  clist_province[[idx]][[years[i-1]]]$data
    for (y in years[(i - to_average):(i-2)])
      mean_previous_provinces = clist_province[[idx]][[y]]$data + mean_previous_provinces
    mean_previous_provinces = mean_previous_provinces/to_average
    mean_previous_provinces = fdata(mean_previous_provinces,argvals)
    basis = create.bspline.basis(rangeval = mean_previous_provinces$rangeval,breaks = knots_vector)
    model = fregre.basis.fr(mean_previous_provinces,clist_province[[idx]][[year]],basis.s = basis,basis.t = basis,lambda.s = 1e-1,lambda.t = 1e-1)
    regr_models_provinces[[year]] = list()
    regr_models_provinces[[year]]$model = model
    regr_models_provinces[[year]]$previous_mean = mean_previous_provinces
    
  }
  
  predictions_province = errors_province = list()
  
  for(i in (to_average+1):(length(years)-1))
  {
    year = years[i+1]
    
    predictions_province[[year]] = predict(regr_models_provinces[[years[i]]]$model, 
                                           regr_models_provinces[[year]]$previous_mean)
    errors_province[[year]] = clist_province[[idx]][[year]] - predictions_province[[year]]
    
    resid = regr_models_provinces[[year]]$model$residuals$data
    y = regr_models_provinces[[year]]$model$y$data
    
    #r.squared for the model given alpha and year
    r.sq[[year]] = 1 - mean(int.simpson(resid^2))/mean(int.simpson((y - mean(y))^2))
    
    meape[[year]] = apply(abs(errors_province[[year]]$data)/
                            abs(clist_province[[idx]][[year]]$data), 1, median)
    
    if (transformation == "Tsagris")
      predictions_province[[year]] = alpha.folding(predictions_province[[year]], alpha)
    
    if (transformation == "Isometric")
      predictions_province[[year]] = alpha.isometric.inv(predictions_province[[year]], alpha)
    
    else 
      predictions_province[[year]] = inv_clr(predictions_province[[year]])
    
  } 
  
  return(list(models = regr_models_provinces, R2 = r.sq, MedAPE = meape, 
              preds = predictions_province, errors = errors_province))
  
}


##### MODELS
load("Output/Data/Smoothing_Provinces_52_alpha_23.Rdata")

alpha.best = 0.85

ncores <- max(parallel::detectCores() -1,1)
if (ncores==1) {
  foreach::registerDoSEQ()
} else{
  cl <- suppressWarnings(parallel::makePSOCKcluster(ncores ))
  doParallel::registerDoParallel(cl)
}
ops.fda.usc(reset=TRUE, verbose=TRUE)

model.tsagris = model.fitting(alpha.best, "Tsagris")
model.clr = model.fitting(0, "CLR")

errors_province_best = model.tsagris$errors
regr_models_provinces_best = model.tsagris$models
predictions_province_best = model.tsagris$preds
r_squared_best = model.tsagris$R2
meape_best = model.tsagris$MedAPE

errors_province_clr = model.clr$errors
regr_models_provinces_clr = model.clr$models
predictions_province_clr = model.clr$preds
r_squared_clr = model.clr$R2
meape_clr = model.clr$MedAPE


save(meape_clr, r_squared_clr, regr_models_provinces_clr, predictions_province_clr, errors_province_clr, file = Rdata_path("Regression_Province_clr_23"))
save(meape_best, r_squared_best, regr_models_provinces_best, predictions_province_best, errors_province_best, alpha.best, file = Rdata_path("Regression_Province_best_alpha_23"))


#############
rm(clist_province)
rm(cl)
closeAllConnections()
load("Output/Data/Smoothing_Provinces_52_alpha_isometric_23.Rdata")

ncores <- max(parallel::detectCores() -1,1)
if (ncores==1) {
  foreach::registerDoSEQ()
} else{
  cl <- suppressWarnings(parallel::makePSOCKcluster(ncores ))
  doParallel::registerDoParallel(cl)
}
ops.fda.usc(reset=TRUE, verbose=TRUE)

alpha.best_isom = 1
model.isometric = model.fitting(alpha.best_isom, "Isometric")

errors_province_best_isom = model.isometric$errors
regr_models_provinces_best_isom = model.isometric$models
predictions_province_best_isom = model.isometric$preds
r_squared_best_isom = model.isometric$R2
meape_best_isom = model.isometric$MedAPE

save(meape_best_isom, r_squared_best_isom, regr_models_provinces_best_isom, predictions_province_best_isom, errors_province_best_isom, alpha.best_isom, file = Rdata_path("Regression_Province_best_alpha_isometric_23"))


######## PLOTS

rm(list=ls())
load("Output/Data/Regression_Province_clr_23.Rdata")
load("Output/Data/Regression_Province_best_alpha_23.Rdata")
load("Output/Data/Regression_Province_best_alpha_isometric_23.Rdata")
load("Output/Data/Smoothing_Provinces_52_alpha_23.Rdata")
source("Scripts/Utility_Functions.R")

to_average = 4
years = paste0("T_",11:23)
years_names = paste0(20,11:23)

comparison = data.frame(
  Year = rep(years_names[(to_average+2):length(years)], each = 3), 
  Transformation = rep(c("CLR","Tsagris","Isometric"), times = length((to_average+2):length(years))),
  R2 = NA 
)

##R2 comparison table
comparison[comparison$Transformation == "CLR",3] = unlist(r_squared_clr)
comparison[comparison$Transformation == "Tsagris",3] = unlist(r_squared_best)
comparison[comparison$Transformation == "Isometric",3] = unlist(r_squared_best_isom)


####### Density comparison plots for KLdiv and MedAPE

constant_density = fdata(rep(1/diff(clist_province[["0"]][["T_11"]]$rangeval), 
                             length(clist_province[["0"]][["T_11"]]$argvals)),
                         clist_province[["0"]][["T_11"]]$argvals)

KL_div = function(errors, baseline_density, transformation, alpha = 0) {
  errors_density = NULL
  
  if (transformation == "CLR") errors_density = inv_clr(errors)
  if (transformation == "Tsagris") errors_density = alpha.folding(errors, alpha)
  if (transformation == "Isometric") errors_density = alpha.isometric.inv(errors, alpha)
  
  return( int.simpson(errors_density * log(errors_density/baseline_density)))
  
}

years_plot = "T_23"
years_plot_name = "2023"
comparison = data.frame(Year = rep(years_plot_name, each = 107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), each=107),
                        KLdiv = NA)
    
comparison[comparison$Transformation == "CLR" & comparison$Year == years_plot_name,3] = 
  KL_div(errors_province_clr[[years_plot]], constant_density, "CLR")

comparison[comparison$Transformation == "$A_{0.85}$" & comparison$Year == years_plot_name,3] =    
  KL_div(errors_province_best[[years_plot]], constant_density, "Tsagris", alpha.best)

comparison[comparison$Transformation == "$A_{1-IT}$" & comparison$Year == years_plot_name,3] = 
  KL_div(errors_province_best_isom[[years_plot]], constant_density, "Isometric", alpha.best_isom)
    
  
g1 = ggplot(comparison %>% filter(Year == "2023"), aes(KLdiv, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("")  +
  theme(legend.position = "none") +
  theme(axis.text.y = element_blank()) 


comparison = data.frame(Year = rep(years_plot_name, each = 107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), each=107), 
                        MedAPE = NA 
)

comparison[comparison$Transformation == "CLR",3] = meape_clr[[years_plot]]
comparison[comparison$Transformation == "$A_{0.85}$",3] = meape_best[[years_plot]]
comparison[comparison$Transformation == "$A_{1-IT}$",3] = meape_best_isom[[years_plot]]

q1 = ggplot(comparison %>% filter(Year == "2023"), aes(MedAPE, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("")  +
  theme(legend.position = "none") +
  theme(axis.text.y = element_blank())


dev.new(width=15, height=7)
(g1 + q1) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_blank())





####### Density comparison plots for KLdiv and MedAPE - SUPPLEMENTARY MATERIAL

constant_density = fdata(rep(1/diff(clist_province[["0"]][["T_11"]]$rangeval), 
                             length(clist_province[["0"]][["T_11"]]$argvals)),
                         clist_province[["0"]][["T_11"]]$argvals)

KL_div = function(errors, baseline_density, transformation, alpha = 0) {
  errors_density = NULL
  
  if (transformation == "CLR") errors_density = inv_clr(errors)
  if (transformation == "Tsagris") errors_density = alpha.folding(errors, alpha)
  if (transformation == "Isometric") errors_density = alpha.isometric.inv(errors, alpha)
  
  return( int.simpson(errors_density * log(errors_density/baseline_density)))
  
}

years_plot = c("T_21", "T_22")
years_plot_name = c("2021", "2022")
comparison = expand.grid(Year = rep(years_plot_name, each = 107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), each=107),
                        KLdiv = NA)

for (i in 1:2) {
    
  comparison[comparison$Transformation == "CLR" & comparison$Year == years_plot_name[i],3] = 
    KL_div(errors_province_clr[[years_plot[i]]], constant_density, "CLR")
  
  comparison[comparison$Transformation == "$A_{0.85}$" & comparison$Year == years_plot_name[i],3] =    
    KL_div(errors_province_best[[years_plot[i]]], constant_density, "Tsagris", alpha.best)
  
  comparison[comparison$Transformation == "$A_{1-IT}$" & comparison$Year == years_plot_name[i],3] = 
    KL_div(errors_province_best_isom[[years_plot[i]]], constant_density, "Isometric", alpha.best_isom)

}

g1 = ggplot(comparison %>% filter(Year == "2021"), aes(KLdiv, colour = Transformation)) +
  geom_density() + xlim(c(0, 0.1524))+ 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2021") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank()) 

g2 = ggplot(comparison %>% filter(Year == "2022"), aes(KLdiv, colour = Transformation)) +
  geom_density() + xlim(c(0, 0.1524))+ 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2022") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank()) 

kldiv = (g1 + g2) + plot_layout(ncol = 2, guides = 'collect') 



comparison = expand.grid(Year = rep(years_plot_name, each = 107),
                         Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), each=107),
                         MedAPE = NA)

comparison[comparison$Transformation == "CLR",3] = c(meape_clr[["T_21"]], meape_clr[["T_22"]])
comparison[comparison$Transformation == "$A_{0.85}$",3] = c(meape_best[["T_21"]], meape_best[["T_22"]])
comparison[comparison$Transformation == "$A_{1-IT}$",3] = c(meape_best_isom[["T_21"]], meape_best_isom[["T_22"]])

q1 = ggplot(comparison %>% filter(Year == "2021"), aes(MedAPE, colour = Transformation)) +
  geom_density() + xlim(c(0.32, 8.66)) +
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2021") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank()) 

q2 = ggplot(comparison %>% filter(Year == "2022"), aes(MedAPE, colour = Transformation)) +
  geom_density() + xlim(c(0.32, 8.66)) +
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2022") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank()) 

medape = (q1 + q2) + plot_layout(ncol = 2, guides = 'collect') 

dev.new(width=15, height=7)
kldiv / (medape & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_blank()))






# 
# ################# confronto medape su stessa scala
# 
# load("Output/Data/Smoothing_Provinces_52_alpha_23.Rdata")
# clist_province_tsag = clist_province
# rm(clist_province)
# 
# load("Output/Data/Smoothing_Provinces_52_alpha_isometric_23.Rdata")
# clist_province_isom = clist_province
# rm(clist_province)
# 
# #### scala CLR
# 
# meape_best_rescaled = meape_best_isom_rescaled = list()
# for (year in names(meape_clr)) {
#   meape_best_rescaled[[year]] = apply(abs(model.tsagris$errors[[year]]$data)/abs(clist_province_tsag[["0"]][[year]]$data),1,median)
#   meape_best_isom_rescaled[[year]] = apply(abs(model.isometric$errors[[year]]$data)/abs(clist_province_tsag[["0"]][[year]]$data),1,median)
# }
# 
# comparison = data.frame(
#   Year = rep(years_names[(to_average+2):length(years)], each = 3 * 107), 
#   Transformation = rep(c("CLR","Tsagris","Isometric"), times = length((to_average+2):length(years)), each=107),
#   MedAPE = NA 
# )
# 
# comparison[comparison$Transformation == "CLR",3] = unlist(meape_clr)
# comparison[comparison$Transformation == "Tsagris",3] = unlist(meape_best_rescaled)
# comparison[comparison$Transformation == "Isometric",3] = unlist(meape_best_isom_rescaled)
# 
# ggplot(comparison, aes(x = as.factor(Year), y = MedAPE, fill = Transformation)) +
#   geom_boxplot(position = position_dodge(width = 0.75)) +  # Align the boxplots
#   geom_vline(xintercept = seq(1.5, length(unique(comparison$Year)) - 0.5, by = 1), linetype = "dashed", color = "grey") +
#   labs(x = "Year", y = "CLR scale", 
#        title = TeX(paste0("MedAPE indicator over the years. Tsagris $\\alpha = ",alpha.best,"$, Isometric $\\alpha = ",alpha.best_isom,"$"))) +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   scale_fill_manual(values = c("CLR" = "skyblue", "Tsagris" = "orange", "Isometric" = "red"))
# 
# 
# 
# #### scala Tsagris
# 
# meape_best_rescaled = meape_best_isom_rescaled = meape_clr_rescaled = list()
# for (year in names(meape_clr)) {
#   meape_clr_rescaled[[year]] = apply(abs(model.clr$errors[[year]]$data)/abs(clist_province_tsag[[paste0(alpha.best)]][[year]]$data),1,median)
#   meape_best_rescaled[[year]] = meape_best[[year]]
#   meape_best_isom_rescaled[[year]] = apply(abs(model.isometric$errors[[year]]$data)/abs(clist_province_tsag[[paste0(alpha.best)]][[year]]$data),1,median)
# }
# 
# comparison = data.frame(
#   Year = rep(years_names[(to_average+2):length(years)], each = 3 * 107), 
#   Transformation = rep(c("CLR","Tsagris","Isometric"), times = length((to_average+2):length(years)), each=107),
#   MedAPE = NA 
# )
# 
# comparison[comparison$Transformation == "CLR",3] = unlist(meape_clr_rescaled)
# comparison[comparison$Transformation == "Tsagris",3] = unlist(meape_best_rescaled)
# comparison[comparison$Transformation == "Isometric",3] = unlist(meape_best_isom_rescaled)
# 
# ggplot(comparison, aes(x = as.factor(Year), y = MedAPE, fill = Transformation)) +
#   geom_boxplot(position = position_dodge(width = 0.75)) +  # Align the boxplots
#   geom_vline(xintercept = seq(1.5, length(unique(comparison$Year)) - 0.5, by = 1), linetype = "dashed", color = "grey") +
#   labs(x = "Year", y = "Tsagris scale", 
#        title = TeX(paste0("MedAPE indicator over the years. Tsagris $\\alpha = ",alpha.best,"$, Isometric $\\alpha = ",alpha.best_isom,"$"))) +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   scale_fill_manual(values = c("CLR" = "skyblue", "Tsagris" = "orange", "Isometric" = "red"))
# 
# 
# 
# #### scala Isometric
# 
# meape_best_rescaled = meape_best_isom_rescaled = meape_clr_rescaled = list()
# for (year in names(meape_clr)) {
#   meape_clr_rescaled[[year]] = apply(abs(model.clr$errors[[year]]$data)/abs(clist_province_isom[[paste0(alpha.best_isom)]][[year]]$data),1,median)
#   meape_best_rescaled[[year]] = apply(abs(model.tsagris$errors[[year]]$data)/abs(clist_province_isom[[paste0(alpha.best_isom)]][[year]]$data),1,median)
#   meape_best_isom_rescaled[[year]] = meape_best_isom[[year]]
# }
# 
# comparison = data.frame(
#   Year = rep(years_names[(to_average+2):length(years)], each = 3 * 107), 
#   Transformation = rep(c("CLR","Tsagris","Isometric"), times = length((to_average+2):length(years)), each=107),
#   MedAPE = NA 
# )
# 
# comparison[comparison$Transformation == "CLR",3] = unlist(meape_clr_rescaled)
# comparison[comparison$Transformation == "Tsagris",3] = unlist(meape_best_rescaled)
# comparison[comparison$Transformation == "Isometric",3] = unlist(meape_best_isom_rescaled)
# 
# ggplot(comparison, aes(x = as.factor(Year), y = MedAPE, fill = Transformation)) +
#   geom_boxplot(position = position_dodge(width = 0.75)) +  # Align the boxplots
#   geom_vline(xintercept = seq(1.5, length(unique(comparison$Year)) - 0.5, by = 1), linetype = "dashed", color = "grey") +
#   labs(x = "Year", y = "Isometric scale", 
#        title = TeX(paste0("MedAPE indicator over the years. Tsagris $\\alpha = ",alpha.best,"$, Isometric $\\alpha = ",alpha.best_isom,"$"))) +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   scale_fill_manual(values = c("CLR" = "skyblue", "Tsagris" = "orange", "Isometric" = "red"))
# 
# 
# 
# ## tutte le scale mi riportano al risultato che Isometric > Tsagris > CLR






