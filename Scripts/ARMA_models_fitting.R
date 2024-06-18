rm(list=ls())
source("Scripts/Utility_Functions.R")
load("Output/Data/Aggregated_Provinces_23.Rdata") ### for smoothing on Provinces

indexes = paste0(seq(0,1,by=0.05))
class = "15_21"
years = paste0("T_", c(11:23))
K = 3

original_discrete = matrix(NA, nrow=107, ncol=0)

for (year in years) {
  temp = provinces_aggregated[,-c(4:29)] %>% filter(CL_ETA == class) %>% ungroup() %>% dplyr::select(NOME_PROVINCIA,partial_date_death,one_of(year)) %>% arrange(partial_date_death)
  
  rows = 1:dim(temp)[1]
  rows = rows[which(temp$partial_date_death != "2024-02-29")]

  temp = temp[rows,]
  temp$partial_date_death = paste0("20",substr(year,3,4),
                                   substr(temp$partial_date_death, 5,10))
  
  temp$dummy = temp[[year]]
  temp = temp %>% dplyr::select(-one_of(year))
  temp = temp %>% pivot_wider(names_from = partial_date_death, values_from = dummy)
  temp[is.na(temp)] = 0
  temp = as.matrix(temp[,-1])

  original_discrete = cbind(original_discrete, temp)  
}

original_discrete = original_discrete/rowSums(original_discrete)

######## smoothing functions

alfa = 5e-3
der = 2
npoints = 10000
knots = 52 * 13

dates = as.numeric(as.Date("2011-01-01")):as.numeric(as.Date("2023-12-31"))
bisestili = as.numeric(c(as.Date("2012-02-29"),as.Date("2016-02-29"),as.Date("2020-02-29")))
dates = dates[-which(dates %in% bisestili)]
xcp = as.numeric(dates)
knots_vector = floor(xcp[seq(1,length(xcp), length.out = knots)])
degree = 3
add_node = 1
add_days = 7

xcp_extended = c((xcp[1] - add_days):(xcp[1] - 1),xcp,(xcp[length(xcp)] + 1):(xcp[length(xcp)] + add_days))  #add one week before ad one week after
knots_vector_extended = c(xcp_extended[1],knots_vector,xcp_extended[length(xcp_extended)])  #54 knots

argvals = seq(xcp_extended[1],xcp_extended[length(xcp_extended)],length.out = npoints)
start = xcp[1]-1
stop = xcp[length(xcp)]+1

alpha.seq = as.numeric(indexes)

smoothing.row <- function(alpha, vector) {
  nbasis <- 4:100
  
  gcv <- numeric(length(nbasis))
  for (j in 1:length(nbasis)){
    basis <- create.fourier.basis(range(argvals), nbasis[j], dropind=1)
    functionalPar <- fdPar(fdobj=basis, Lfdobj=der, lambda=alfa)  
    gcv[j] <- smooth.basis(xcp_extended, vector, functionalPar)$gcv
  }
  
  basis <- create.fourier.basis(rangeval=range(argvals), nbasis=nbasis[which.min(gcv)], dropind=1)
  functionalPar <- fdPar(fdobj=basis, Lfdobj=der, lambda=alfa)  
  Xsp <- smooth.basis(argvals=xcp_extended, y=vector, fdParobj=functionalPar)
  ret <- eval.fd(evalarg = argvals, Xsp$fd) #  the curve smoothing the data
  
  t(ret)
}

smoothing <- function(alpha, matrix) {
  
  ret = NULL
  
  if (alpha == 0) {  
    smooth = smoothSplines(k = degree,l = der,alpha = alfa,data = matrix,
                           knots = knots_vector_extended, num_points = npoints,xcp = xcp_extended, fast = 1)
    smoothed = fdata(smooth$Y, argvals = argvals)
    smoothed = conditioning(smoothed, start, stop)
    ret = clr(smoothed)
  }
  
  else {
    
    matrix.alpha = alpha.tsagris(fdata(matrix), alpha)$data
    
    z = matrix(nrow=107, ncol=npoints)
    foreach(i=1:107) %dopar% {
      smoothing.row(alpha, matrix.alpha[i,])
    } -> result
    
    # Combine results into matrix z
    for (i in 1:107) {
      z[i,] <- result[[i]]
    }
    
    z = fdata(z, argvals = argvals)
    indices = which(z$argvals > start & z$argvals< stop)
    ret = fdata(z$data[,indices], z$argvals[indices])
  }
  
  ret
  
}

ncores <- max(parallel::detectCores() -1,1)
cl <- makeCluster(ncores)
clusterExport(cl, 
              c("xcp_extended", "alfa", "der", "degree", "knots_vector_extended",
                "smoothing.row", "argvals", "xcp_extended","create.fourier.basis", 
                "smooth.basis", "eval.fd", "fdPar","npoints"))
registerDoParallel(cl)

############

original_transformed = list()

for (alpha in alpha.seq) {
  
  print("---------------")
  idx = paste0(alpha)
  print(idx)
  
  matrix = cbind(original_discrete[,length(xcp):(length(xcp) - add_days + 1)], 
                 original_discrete, 
                 original_discrete[,1:add_days])
  
  original_transformed[[idx]] = smoothing(alpha, matrix)
  
}

original_density = matrix(0, nrow=107, ncol=ncol(original_transformed[["0"]]$data))
for (alpha in alpha.seq[-1]) {
  to_add = alpha.folding(original_transformed[[paste0(alpha)]], alpha)
  original_density = original_density + to_add$data/20
}
original_density = fdata(original_density, original_transformed[["0"]]$argvals)

save(original_density, original_transformed, file = Rdata_path("density_2023") )

###########

load("Output/Data/density_2023.Rdata")

test_series = time_series_PC = list()
columns_23 = which(as.Date(original_density$argvals, origin = "1970-01-01") > 
                     as.Date("2022-12-31", origin = "1970-01-01" ))

##build a time series considering all the years together for 70+ age class
for (alpha_index in indexes) {
  
  print(alpha_index)
  
  data_matrix = fdata(original_transformed[[alpha_index]]$data[,-columns_23],
                      original_transformed[[alpha_index]]$argvals[-columns_23])

  ##functional PCA
  time_series_PC[[alpha_index]] =
    fdata2pc(data_matrix, ncomp = K, norm = F)
  
  test_series[[alpha_index]] = fdata(original_transformed[[alpha_index]]$data[,columns_23],
                                     original_transformed[[alpha_index]]$argvals[columns_23])
  
}

save(time_series_PC, test_series, file = Rdata_path("PC_Province_alpha_TS"))

########################

load("Output/Data/PC_Province_alpha_TS.Rdata") #without repeating all of the above code

indexes = paste0(seq(0,1,by=0.05))
class = "15_21"  #for simplicity
years = paste0("T_", seq(11,23))

#ARMA model for the means
models_mean = list()
seasonality_mean = list() 

for (alpha_index in indexes) {
  
  print(alpha_index)
  
  mean_func = time_series_PC[[alpha_index]]$mean$data[1,-1]
  
  ##create time series objects
  ts_func = ts(mean_func, start=c(2011,1), end = c(2022,767), frequency = 767)
  
  seasonality_mean[[alpha_index]] = decompose(ts_func, "additive")$seasonal
  
  ##remove seasonality
  adjusted = ts_func - seasonality_mean[[alpha_index]]
  
  models_mean[[alpha_index]] = auto.arima(x = adjusted, method="CSS", seasonal = F, 
                                          stationary = F, allowdrift=F)
  
}


#ARMA models for the first K eigenfunctions
models_eigen = list()
seasonality_eigen = list() 

for (alpha_index in indexes) {
  
  print("-----------------")
  print(alpha_index)
  
  models_eigen[[alpha_index]] = list()
  seasonality_eigen[[alpha_index]] = list() 
  
  eigenfunctions = time_series_PC[[alpha_index]]$rotation$data[1:K,]
  
  for (i in 1:K) {
    
    print(i)
    
    ##create time series objects
    ts_func = ts(eigenfunctions[i,], start=c(2011,1), end = c(2022,767), frequency = 767)
    
    seasonality_eigen[[alpha_index]][[i]] = decompose(ts_func, "additive")$seasonal
    
    ##remove seasonality
    adjusted = ts_func - seasonality_eigen[[alpha_index]][[i]]
    
    models_eigen[[alpha_index]][[i]] = 
      auto.arima(x = adjusted, method="CSS", seasonal = F, stationary = F, allowdrift=F)
    
  }
  
  
}

rm(ts_func, adjusted, eigenfunctions, mean_func)

save(models_mean, models_eigen, seasonality_eigen, seasonality_mean, file = Rdata_path("ARMA_PCA_models_alpha"))


######################## PERFORMANCE EVALUATION

rm(list=ls())
load("Output/Data/ARMA_PCA_models_alpha.Rdata")
load("Output/Data/PC_Province_alpha_TS.Rdata") 
load("Output/Data/Aggregated_Provinces_23.Rdata")
load("Output/Data/density_2023.Rdata")
source("Scripts/Utility_Functions.R")

indexes = paste0(seq(0,1,by=0.05))
class = "15_21"  #for simplicity
years = paste0("T_", seq(11,23))
K = 3

stop_idx = max(which(as.Date(original_density$argvals, origin="1970-01-01")<as.Date("2023-01-01")))+1
stop = original_density$argvals[stop_idx]+1

#from 2022/01/01
# start_idx = min(which(as.Date(original_density$argvals, origin="1970-01-01")>as.Date("2022-01-01")))-1
# start = original_density$argvals[start_idx]

# #from 2011/01/01
start_idx = min(which(as.Date(original_density$argvals, origin="1970-01-01")>as.Date("2011-01-01")))-1
start = original_density$argvals[start_idx]

##conditional density 
conditional_density = conditioning(original_density, start, stop)

##MedAPE, KL-divergence and RRSE in L2
meape = rrse = divergence = list()

##mean of KL-divergence for the selection of the best alpha
indicator = list()

for (alpha_index in indexes) {
  
  #total seasonality effect
  seasonality = matrix( rep(seasonality_mean[[alpha_index]][768:(767*2)], 107), ncol = 767, byrow = T)
  
  #compute predictions in each province 
  preds = matrix(NA, 107, 767)
  scores = time_series_PC[[alpha_index]]$x[,1:K]
  
  for (i in 1:107) {
    
    preds[i,] = predict(models_mean[[alpha_index]], n.ahead = 767, se.fit=F)
    
    for (j in 1:K) {

      preds[i,] = preds[i,] + scores[i,j] *
        predict(models_eigen[[alpha_index]][[j]], n.ahead=767, se.fit=F)

      seasonality[i,] = seasonality[i,] + scores[i,j] * seasonality_eigen[[alpha_index]][[j]][768:(767*2)]

    }

  }
  
  alpha = as.numeric(alpha_index)
  
  ##function in L2 made of original data up to 2022 
  predicted_transformed = original_transformed[[alpha_index]]$data[,1:(stop_idx-1)]
  
  ##add predictions for 2023 and convert to fdata object
  predicted_transformed = cbind(predicted_transformed, preds+seasonality)
  predicted_transformed = fdata(predicted_transformed, original_density$argvals)
  
  ##resulting density
  predicted_density = alpha.folding(predicted_transformed, alpha)
  
  ##resulting conditional density
  predicted_conditional_density = conditioning(predicted_density, start, stop)
  
  true = alpha.tsagris(conditional_density,alpha)$data
  pred = alpha.tsagris(predicted_conditional_density, alpha)$data
  
  meape[[alpha_index]] = rrse[[alpha_index]] = numeric(107)
  
  divergence[[alpha_index]] = int.simpson(predicted_conditional_density * log(predicted_conditional_density/conditional_density))
  
  for (k in 1:107) {
    meape[[alpha_index]][k] = MedianAPE(y_pred = pred[k,], y_true = true[k,])
    rrse[[alpha_index]][k] = RRSE(y_pred = pred[k,], y_true = true[k,])
  }
  
  indicator[[alpha_index]] = mean(divergence[[alpha_index]])
}


########### plots
data = expand.grid(Alpha = rep(indexes, each=107), 
                   KLdiv = NA,
                   MedAPE = NA,
                   RRSE = NA)

for (idx in indexes) {
  data$KLdiv[data$Alpha == idx] = divergence[[idx]]
  data$MedAPE[data$Alpha == idx] = meape[[idx]]
  data$RRSE[data$Alpha == idx] = rrse[[idx]]
}

kl = ggplot(data, aes(x = as.factor(Alpha), y = KLdiv)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape =NA, fill="lightblue") +  
  labs(y = "", x = TeX("$\\alpha$")) +
  theme_minimal() + ggtitle("KL divergence") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) +
  ylim(c(0,0.0125)) #long time span 
  #ylim(c(0,0.01)) #short time span

ma = ggplot(data, aes(x = as.factor(Alpha), y = MedAPE)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape =NA, fill="lightgreen") +  
  labs(y = "", x = TeX("$\\alpha$")) +
  theme_minimal() + ggtitle("MedAPE") + 
  theme(legend.position = "none" , plot.title = element_text(hjust = 0.5)) +
  ylim(c(0,1.25)) #both time spans

rr = ggplot(data, aes(x = as.factor(Alpha), y = RRSE)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape =NA, fill="orange") +  
  labs(y = "", x = TeX("$\\alpha$")) +
  theme_minimal() + ggtitle("RRSE") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) + 
  ylim(c(0,1.25)) #both time spans

dev.new(width=12, height=7.5)
((kl | ma) / rr)  #boxplots of the indicators depending on the value of alpha


idx.best = indexes[which.min(indicator)]
alpha.best = as.numeric(idx.best)

######

best_model = list(mean = models_mean[[idx.best]],
                  eigen = models_eigen[[idx.best]],
                  alpha.best = as.numeric(idx.best),
                  meape = meape[[idx.best]],
                  rrse = rrse[[idx.best]])

clr_model = list(mean = models_mean[["0"]],
                 eigen = models_eigen[["0"]],
                 meape = meape[["0"]],
                 rrse = rrse[["0"]])

save(best_model, clr_model, file=Rdata_path("ARMA_models_23_best"))

closeAllConnections()




## 0.4 è l'alpha ottimo per le condizionali 2011-2022 
## 0.3 è l'alpha ottimo per le condizionali 2023 ma non ci sono evidenti differenze 


