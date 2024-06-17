rm(list=ls())
source("Scripts/Utility_Functions.R")

# Regression models on Provinces
# Create the regression models for each year and class. Set some parameters and utilities:

classe = c("0_10","11_14","15_21")
years = paste0("T_",11:20)
knots = 52
xcp = as.numeric((as.numeric(as.Date("2020-01-01")):as.numeric(as.Date("2020-12-31"))))
knots_vector = floor(xcp[seq(1,length(xcp), length.out = knots)])
to_average = 4
alpha.seq = seq(0, 1, by=0.05)  ##alpha transformation parameter
class = classe[3]

#######FUNCTIONS

models.fitting = function(alpha, train.indexes, test.indexes) {
  idx = paste(alpha,"", sep="")
  
  regr_models_provinces = r.sq = meape =  list()
  
  #fit the models
  for(i in (to_average+1):length(years))
  {
    year = years[i]
    print(year)
    argvals = clist_province[[idx]][[year]]$argvals
    
    mean_previous_provinces =  clist_province[[idx]][[years[i-1]]]$data
    for (y in years[(i - to_average):(i-2)]) 
      mean_previous_provinces = clist_province[[idx]][[y]]$data + mean_previous_provinces
    
    mean_previous_provinces = mean_previous_provinces/to_average
    
    test.regressor = fdata(mean_previous_provinces[test.indexes,], argvals)
    mean_previous_provinces = fdata(mean_previous_provinces[train.indexes,],argvals)
    
    basis = create.bspline.basis(rangeval = mean_previous_provinces$rangeval,breaks = knots_vector)
    training.set = fdata(clist_province[[idx]][[year]]$data[train.indexes,], argvals)
    
    model = fregre.basis.fr(mean_previous_provinces, training.set, basis.s = basis, basis.t = basis,lambda.s = 1e-1,lambda.t = 1e-1)
    regr_models_provinces[[year]] = list()
    regr_models_provinces[[year]]$model = model
    regr_models_provinces[[year]]$previous_mean = mean_previous_provinces
    
    #compute the r.squared coefficient for each year 
    test.set = clist_province[[idx]][[year]]$data[test.indexes,]
    
    r.sq[[year]] = 1 - mean(int.simpson(model$residuals$data^2))/mean(int.simpson((model$y$data - mean(model$y$data))^2))
    regr_models_provinces[[year]]$test.set = test.set
    regr_models_provinces[[year]]$test.regressor = test.regressor
    
  }
  
  return(list(models = regr_models_provinces, R2 = r.sq))
}


cross.validation.indicators = function(nfolds=5, alpha, seed=2024, transformation = "Tsagris") {
  
  set.seed(seed)
  idx = paste0(alpha)
  print(idx)
  
  rows = sample(1:107, 107, replace=F)
  folds = suppressWarnings(split(rows,1:nfolds))
  
  ##R2 for each year: one value per each fold
  ##MedAPE for each year: one value for each fold = mean of MedAPEs of out of bag province 
  meape = r.squared = errors_divergence = list()
  
  foreach(j=1:nfolds, .packages = "fda.usc", .inorder = T) %dopar% {
    test.indexes = folds[[j]]
    train.indexes = (1:107)[-test.indexes]
    
    models.fitting(alpha, train.indexes, test.indexes)
    
  }  -> cv.models
  
  for (y in names(cv.models[[1]]$R2)) {
    
    r.squared[[y]] = meape[[y]] = errors_divergence[[y]] = numeric(nfolds)
    
    for (k in 1:nfolds) {
      r.squared[[y]][[k]] = cv.models[[k]]$R2[[y]]
      model = cv.models[[k]]$models[[y]]$model
      test.regressor = cv.models[[k]]$models[[y]]$test.regressor
      test.set = cv.models[[k]]$models[[y]]$test.set
      predictions = predict(model, test.regressor)$data
      errors = test.set - predictions
      meape[[y]][[k]] = mean(apply(abs(errors)/abs(test.set), 1, median))
      
      errors_density = NULL
      
      if (transformation == "Tsagris") errors_density = alpha.folding(fdata(errors, constant_density$argvals),alpha)
      else errors_density = alpha.isometric.inv(fdata(errors, constant_density$argvals),alpha)
      
      errors_divergence[[y]][[k]] = mean(int.simpson(errors_density * log(errors_density/constant_density)))
    }
  }
  
  return(list(R2 = r.squared, MedAPE = meape, KLdiv = errors_divergence))
  
}


####### TSAGRIS

load("Output/Data/Smoothing_Provinces_52_alpha.Rdata")

constant_density = fdata(rep(1/diff(clist_province[["0"]][["T_11"]]$rangeval), 
                             length(clist_province[["0"]][["T_11"]]$argvals)),
                         clist_province[["0"]][["T_11"]]$argvals)

nfolds = 5
indicators_data = list()
years_names = paste0(20, substr(years,3,4))

ncores <- max(parallel::detectCores() -1,1)
if (ncores==1) {
  foreach::registerDoSEQ()
} else{
  cl <- suppressWarnings(parallel::makePSOCKcluster(ncores ))
  clusterExport(cl, list("knots_vector","models.fitting","to_average","years","clist_province",
                         "create.bspline.basis"))
  registerDoParallel(cl)
}
ops.fda.usc(reset=TRUE, verbose=TRUE)

for (alpha in alpha.seq) {
  
  idx = paste0(alpha)
  data = data.frame(Year = rep(years_names[(to_average+2):length(years)], each=nfolds*3), 
                    Indicator = rep(c("R^2", "MedAPE","KL-divergence"), each=nfolds, times=length(years)-to_average-1),
                    Value = NA)
  
  cv.indicators = cross.validation.indicators(nfolds, alpha)
  
  for (year in years) {
    name = paste0(20, substr(year,3,4))
    data[data$Year == name & data$Indicator == "R^2",3] = cv.indicators$R2[[year]]
    data[data$Year == name & data$Indicator == "MedAPE",3] = cv.indicators$MedAPE[[year]]
    data[data$Year == name & data$Indicator == "KL-divergence",3] = cv.indicators$KLdiv[[year]]
  }
  
  indicators_data[[idx]] = data
  
}


save(indicators_data, file = Rdata_path("Indicators_cv_Tsagris_20"))


load("Output/Data/Indicators_cv_Tsagris_20.Rdata")

kl_div = NULL

for (idx in paste0(alpha.seq)) {
  
  kl_div_mean = indicators_data[[idx]] %>% filter(Indicator == "KL-divergence") %>% group_by(Year) %>%
    summarise(mean = mean(Value))
  
  kl_div = c(kl_div, mean(kl_div_mean$mean))
  
}

names(kl_div) = alpha.seq
idx.best = names(kl_div)[which.min(kl_div)] #0.85

closeAllConnections()



####### ISOMETRIC
rm(clist_province)
load("Output/Data/Smoothing_Provinces_52_alpha_isometric.Rdata")

constant_density = fdata(rep(1/diff(clist_province[["0"]][["T_11"]]$rangeval), 
                             length(clist_province[["0"]][["T_11"]]$argvals)),
                         clist_province[["0"]][["T_11"]]$argvals)
nfolds = 5
indicators_data_isom = list()
years_names = paste0(20, substr(years,3,4))

ncores <- max(parallel::detectCores()-1,1)
if (ncores==1) {
  foreach::registerDoSEQ()
} else{
  cl <- suppressWarnings(parallel::makePSOCKcluster(ncores ))
  clusterExport(cl, list("knots_vector","models.fitting","to_average","years","clist_province",
                         "create.bspline.basis"))
  registerDoParallel(cl)
}
ops.fda.usc(reset=TRUE, verbose=TRUE)

for (alpha in alpha.seq) {
  
  idx = paste0(alpha)
  data = data.frame(Year = rep(years_names[(to_average+2):length(years)], each=nfolds*3), 
                    Indicator = rep(c("R^2", "MedAPE","KL-divergence"), each=nfolds, times=length(years)-to_average-1),
                    Value = NA)
  
  cv.indicators = cross.validation.indicators(nfolds = nfolds, alpha = alpha, transformation = "Isometric")
  
  for (year in years) {
    name = paste0(20, substr(year,3,4))
    data[data$Year == name & data$Indicator == "R^2",3] = cv.indicators$R2[[year]]
    data[data$Year == name & data$Indicator == "MedAPE",3] = cv.indicators$MedAPE[[year]]
    data[data$Year == name & data$Indicator == "KL-divergence",3] = cv.indicators$KLdiv[[year]]
  }
  
  indicators_data_isom[[idx]] = data
  
}

save(indicators_data_isom, file = Rdata_path("Indicators_cv_Isometric_20"))


load("Output/Data/Indicators_cv_Isometric_20.Rdata")


kl_div = NULL

for (idx in paste0(alpha.seq)) {
  
  kl_div_mean = indicators_data_isom[[idx]] %>% filter(Indicator == "KL-divergence") %>% group_by(Year) %>%
    summarise(mean = mean(Value))
  
  kl_div = c(kl_div, mean(kl_div_mean$mean))
  
}

names(kl_div) = alpha.seq
idx.best_isom = names(kl_div)[which.min(kl_div)]  #1


############################  comparison between CLR, Tsagris and Isometric in CV performance

rm(list=ls())
load("Output/Data/Indicators_cv_Tsagris_20.Rdata")
load("Output/Data/Indicators_cv_Isometric_20.Rdata")

comparison = rbind(indicators_data[["0"]], 
                   indicators_data[["0.85"]], 
                   indicators_data_isom[["1"]])
comparison$Alpha = c(rep(0,75),rep(0.85,75),rep(1,75))
comparison$Transformation = c(rep("CLR",75), rep("$A_{0.85}$",75), rep("$A_{1-IT}$",75) )

temp = comparison[comparison$Indicator == "KL-divergence",]
temp$Transformation = as.factor(temp$Transformation)

dev.new(width=10, height=5)
ggplot(temp, aes(x = as.factor(Year), y = Value, fill = Transformation)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Align the boxplots
  geom_vline(xintercept = seq(1.5, length(unique(temp$Year)) - 0.5, by = 1), linetype = "dashed", color = "grey") +
  labs(x = "Year", y = "KL-divergence") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
  scale_fill_manual(values = c("CLR" = "skyblue", "$A_{0.85}$" = "orange", "$A_{1-IT}$" = "red"),
                    labels = c("CLR" = TeX("CLR"), "$A_{0.85}$" = TeX("$A_{0.85}$"), "$A_{1-IT}$" = TeX("$A_{1-IT}$")))

