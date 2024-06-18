rm(list=ls())
source("Scripts/Utility_Functions.R")

classe = c("0_10","11_14","15_21")
years = paste0("T_",11:20)
years_names = paste0(20,11:20)
knots = 52
xcp = as.numeric((as.numeric(as.Date("2020-01-01")):as.numeric(as.Date("2020-12-31"))))
knots_vector = floor(xcp[seq(1,length(xcp), length.out = knots)])
to_average = 4
alpha.seq = seq(0, 1, by=0.05)  ##alpha transformation parameter
class = classe[3]

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
load("Output/Data/Smoothing_Provinces_52_alpha.Rdata")

alpha.best = 0.85

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


save(meape_clr, r_squared_clr, regr_models_provinces_clr, predictions_province_clr, errors_province_clr, file = Rdata_path("Regression_Province_clr"))
save(meape_best, r_squared_best, regr_models_provinces_best, predictions_province_best, errors_province_best, alpha.best, file = Rdata_path("Regression_Province_best_alpha"))


load("Output/Data/Smoothing_Provinces_52_alpha_isometric.Rdata")

alpha.best_isom = 1

model.isometric = model.fitting(alpha.best_isom, "Isometric")

errors_province_best_isom = model.isometric$errors
regr_models_provinces_best_isom = model.isometric$models
predictions_province_best_isom = model.isometric$preds
r_squared_best_isom = model.isometric$R2
meape_best_isom = model.isometric$MedAPE

save(meape_best_isom, r_squared_best_isom, regr_models_provinces_best_isom, predictions_province_best_isom, errors_province_best_isom, alpha.best_isom, file = Rdata_path("Regression_Province_best_alpha_isometric"))


######## PLOTS

rm(list=ls())
load("Output/Data/Regression_Province_clr.Rdata")
load("Output/Data/Regression_Province_best_alpha.Rdata")
load("Output/Data/Regression_Province_best_alpha_isometric.Rdata")
load("Output/Data/Smoothing_Provinces_52_alpha.Rdata")
source("Scripts/Utility_Functions.R")

years = paste0("T_",11:20)
years_names = paste0(20,11:20)
to_average = 4

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

years_plot = c("T_16", "T_18","T_20")
years_plot_name = c("2016", "2018","2020")
comparison = data.frame(Year = rep(years_plot_name, each = 3*107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), times = 3, each=107),
                        KLdiv = NA)

for (i in (to_average+2):length(years)) {
  
  year = years_names[i]
  name = years[i]
  
  if (name %in% years_plot) {
    
    comparison[comparison$Transformation == "CLR" & comparison$Year == year,3] = 
      KL_div(errors_province_clr[[name]], constant_density, "CLR")
    
    comparison[comparison$Transformation == "$A_{0.85}$" & comparison$Year == year,3] =    
      KL_div(errors_province_best[[name]], constant_density, "Tsagris", alpha.best)
  
    comparison[comparison$Transformation == "$A_{1-IT}$" & comparison$Year == year,3] = 
      KL_div(errors_province_best_isom[[name]], constant_density, "Isometric", alpha.best_isom)
    
  }
}

g1 = ggplot(comparison %>% filter(Year == "2016"), aes(KLdiv, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                    labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2016") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank()) 
  
g2 = ggplot(comparison %>% filter(Year == "2018"), aes(KLdiv, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2018") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())

g3 = ggplot(comparison %>% filter(Year == "2020"), aes(KLdiv, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2020") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())

kldiv = (g1 + g2 + g3) + plot_layout(ncol = 3, guides = 'collect') 


comparison = data.frame(
  Year = rep(years_names[(to_average+2):length(years)], each = 3 * 107), 
  Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), times = length((to_average+2):length(years)), each=107),
  MedAPE = NA 
)

comparison[comparison$Transformation == "CLR",3] = unlist(meape_clr)
comparison[comparison$Transformation == "$A_{0.85}$",3] = unlist(meape_best)
comparison[comparison$Transformation == "$A_{1-IT}$",3] = unlist(meape_best_isom)


q1 = ggplot(comparison %>% filter(Year == "2016"), aes(MedAPE, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2016") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())

q2 = ggplot(comparison %>% filter(Year == "2018"), aes(MedAPE, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2018") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())


q3 = ggplot(comparison %>% filter(Year == "2020"), aes(MedAPE, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2020") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())


medape = (q1 + q2 + q3) + plot_layout(ncol = 3, guides = 'collect') 
 
dev.new(width=15, height=7)
kldiv / (medape & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_blank()))


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

years_plot = c("T_17", "T_19")
years_plot_name = c("2017", "2019")
comparison = data.frame(Year = rep(years_plot_name, each = 3*107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), times = 2, each=107),
                        KLdiv = NA)

for (i in (to_average+2):length(years)) {
  
  year = years_names[i]
  name = years[i]
  
  if (name %in% years_plot) {
    
    comparison[comparison$Transformation == "CLR" & comparison$Year == year,3] = 
      KL_div(errors_province_clr[[name]], constant_density, "CLR")
    
    comparison[comparison$Transformation == "$A_{0.85}$" & comparison$Year == year,3] =    
      KL_div(errors_province_best[[name]], constant_density, "Tsagris", alpha.best)
    
    comparison[comparison$Transformation == "$A_{1-IT}$" & comparison$Year == year,3] = 
      KL_div(errors_province_best_isom[[name]], constant_density, "Isometric", alpha.best_isom)
    
  }
}

g1 = ggplot(comparison %>% filter(Year == "2017"), aes(KLdiv, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2017") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank()) 

g2 = ggplot(comparison %>% filter(Year == "2019"), aes(KLdiv, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("KL divergence") + ylab("") + ggtitle("2019") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())

kldiv = (g1 + g2) + plot_layout(ncol = 2, guides = 'collect') 


comparison = data.frame(Year = rep(years_plot_name, each = 3*107),
                        Transformation = rep(c("CLR","$A_{0.85}$","$A_{1-IT}$"), times = 2, each=107),
                        MedAPE = NA)

comparison[comparison$Transformation == "CLR",3] = c(meape_clr$T_17, meape_clr$T_19)
comparison[comparison$Transformation == "$A_{0.85}$",3] = c(meape_best$T_17, meape_best$T_19)
comparison[comparison$Transformation == "$A_{1-IT}$",3] = c(meape_best_isom$T_17, meape_best_isom$T_19)


q1 = ggplot(comparison %>% filter(Year == "2017"), aes(MedAPE, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2017") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())

q2 = ggplot(comparison %>% filter(Year == "2019"), aes(MedAPE, colour = Transformation)) +
  geom_density() + 
  scale_colour_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                      labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$"))) +
  xlab("MedAPE") + ylab("") + ggtitle("2019") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_blank())


medape = (q1 + q2) + plot_layout(ncol = 2, guides = 'collect') 

dev.new(width=15, height=7)
kldiv / (medape & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_blank()))



# 
# ################# Confronto MedAPE su stessa scala
# 
# load("Output/Data/Smoothing_Provinces_52_alpha.Rdata")
# clist_province_tsag = clist_province
# rm(clist_province)
# 
# load("Output/Data/Smoothing_Provinces_52_alpha_isometric.Rdata")
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
# 
# 
# 
# 


