rm(list=ls())
source("Scripts/Utility_Functions.R")
load("Output/Data/Aggregated_Provinces_24.Rdata") 

indexes = paste0(seq(0,1,by=0.05))
class = "15_21"
years = paste0("T_", c(11:24))
K = 3

original_discrete = matrix(NA, nrow=107, ncol=0)

for (year in years) {
  temp = provinces_aggregated[,-c(4:29)] %>% filter(CL_ETA == class) %>% ungroup() %>% 
    dplyr::select(NOME_PROVINCIA,partial_date_death,one_of(year)) %>% arrange(partial_date_death)
  
  rows = 1:dim(temp)[1]
  rows = rows[which(temp$partial_date_death != "2024-02-29")]
  
  if (year == "T_24") rows = rows[which(temp$partial_date_death < "2024-02-29")]
  
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
npoints = 12000
knots = 52 * 13 + 8

dates = as.numeric(as.Date("2011-01-01")):as.numeric(as.Date("2024-02-28"))
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

alpha_best = 0.4

smoothing.row <- function(vector) {
  nbasis <- 4:110
  
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
      smoothing.row(matrix.alpha[i,])
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

####################

original_density = original_transformed = list()

for (alpha in c(0, alpha_best)) {
  idx = paste0(alpha)
  
  matrix = cbind(original_discrete[,length(xcp):(length(xcp) - add_days + 1)], 
                 original_discrete, 
                 original_discrete[,1:add_days])
  
  original_transformed[[idx]] = smoothing(alpha, matrix)
  
  original_density[[idx]] = alpha.folding(original_transformed[[idx]], alpha)
  
}

save(original_density, original_transformed, file = Rdata_path("density_2024") )

###########

load("Output/Data/density_2024.Rdata")

columns_24 = which(as.Date(original_density[["0"]]$argvals) > as.Date("2023-12-31"))

##build a time series considering all the years together for 70+ age class

time_series_PC = test_series = list()

for (alpha in c(0, alpha_best)) {
  idx = paste0(alpha)
  data_matrix = fdata(original_transformed[[idx]]$data[,-columns_24],
                    original_transformed[[idx]]$argvals[-columns_24])

  ##functional PCA
  time_series_PC[[idx]] = fdata2pc(data_matrix, ncomp = K, norm = F)
  
  test_series[[idx]] = fdata(original_transformed[[idx]]$data[,columns_24],
                             original_transformed[[idx]]$argvals[columns_24])
}

save(time_series_PC, test_series, file = Rdata_path("PC_Province_alpha_TS_24"))

#################

load("Output/Data/PC_Province_alpha_TS_24.Rdata")

#ARMA models for the mean
models_mean = seasonality_mean = list()

for (alpha in c(0, alpha_best)) {
  idx = paste0(alpha)
  
  mean_func = time_series_PC[[idx]]$mean$data[1,]
  
  ##create time series objects
  ts_func = ts(mean_func, start=c(2011,1), end = c(2023,910), frequency = 910)
  
  seasonality_mean[[idx]] = decompose(ts_func, "additive")$seasonal
  
  ##remove seasonality
  adjusted = ts_func - seasonality_mean[[idx]]
  
  models_mean[[idx]] = auto.arima(x = adjusted, method="CSS", seasonal = F, stationary = F, allowdrift=F)
  
}



#ARMA models for the first K eigenfunctions
seasonality_eigen = models_eigen = list()

for (alpha in c(0, alpha_best)) {
  
  idx = paste0(alpha)
  eigenfunctions = time_series_PC[[idx]]$rotation$data[1:K,]
  
  seasonality_eigen[[idx]] = models_eigen[[idx]] = list()
  
  for (i in 1:K) {
    
    print(i)
    ##create time series objects
    ts_func = ts(eigenfunctions[i,], start=c(2011,1), end = c(2023,910), frequency = 910)
    
    seasonality_eigen[[idx]][[i]] = decompose(ts_func, "additive")$seasonal
    
    ##remove seasonality
    adjusted = ts_func - seasonality_eigen[[idx]][[i]]
    
    models_eigen[[idx]][[i]] = 
      auto.arima(x = adjusted, method="CSS", seasonal = F, stationary = F, allowdrift=F)
    
  }
}


save(models_mean, models_eigen, seasonality_eigen, seasonality_mean, file = Rdata_path("ARMA_PCA_models_24"))


#################### 

rm(list=ls())
load("Output/Data/PC_Province_alpha_TS_24.Rdata") #without repeating all of the above code
load("Output/Data/ARMA_PCA_models_24.Rdata") 
load("Output/Data/ARMA_models_23_best.Rdata")  #load the best model of 2023
load("Output/Data/Aggregated_Provinces_24.Rdata")
load("Output/Data/density_2024.Rdata")
source("Scripts/Utility_Functions.R")

indexes = paste0(seq(0,1,by=0.05))
class = "15_21"  #for simplicity
years = paste0("T_", seq(11,24))
K = 3
alpha.best = 0.4

stop_idx = max(which(as.Date(original_density[["0"]]$argvals, origin="1970-01-01")<as.Date("2024-03-01")))
stop = original_density[["0"]]$argvals[stop_idx]+1

#from 2024/01/01
start_idx = max(which(as.Date(original_density[["0"]]$argvals, origin="1970-01-01")<as.Date("2024-01-01")))
start = original_density[["0"]]$argvals[start_idx]


###### performance evaluation

##MedAPE, KL-divergence and RRSE in L2
meape = rrse = divergence = list()

#save the forecasting residuals to plot them 
residuals = list()

for (alpha in c(0, alpha.best)) {
  idx = paste0(alpha)
  
  ##conditional density
  conditional_density = conditioning(original_density[[idx]], start, stop)

  #total seasonality effect
  seasonality = matrix( rep(seasonality_mean[[idx]][911:(910*2)], 107), ncol = 910, byrow = T)
  
  #compute predictions in each province 
  preds = matrix(NA, 107, 910)
  scores = time_series_PC[[idx]]$x[,1:K]
  
  for (i in 1:107) {
    
    preds[i,] = predict(models_mean[[idx]], n.ahead = 910, se.fit=F)
    
    for (j in 1:K) {
      
      preds[i,] = preds[i,] + scores[i,j] *
        predict(models_eigen[[idx]][[j]], n.ahead=910, se.fit=F)
      
      seasonality[i,] = seasonality[i,] + scores[i,j] * seasonality_eigen[[idx]][[j]][911:(910*2)]
      
    }
    
  }
  
  window = 1:length(conditional_density$argvals)
  
  pred = fdata((preds+seasonality)[,window], conditional_density$argvals) #predictions in L2

  ##resulting conditional density
  predicted_conditional_density = alpha.folding(pred, alpha)
  
  true = alpha.tsagris(conditional_density,alpha)
  
  residuals[[idx]] = true-pred 

  meape[[idx]] = rrse[[idx]] = numeric(107)
  
  divergence[[idx]] = int.simpson(predicted_conditional_density * log(predicted_conditional_density/conditional_density))
  
  for (k in 1:107) {
    meape[[idx]][k] = MedianAPE(y_pred = pred$data[k,], y_true = true$data[k,])
    rrse[[idx]][k] = RRSE(y_pred = pred$data[k,], y_true = true$data[k,])
  }
  
  
}



########### boxplots of predictive performance

data = expand.grid(Alpha = rep(paste0(c(0,alpha.best)), each=107), 
                   KLdiv = NA,
                   MedAPE = NA,
                   RRSE = NA)

for (idx in paste0(c(0,alpha.best))) {
  data$KLdiv[data$Alpha == idx] = divergence[[idx]]
  data$MedAPE[data$Alpha == idx] = meape[[idx]]
  data$RRSE[data$Alpha == idx] = rrse[[idx]]
}

kl = ggplot(data, aes(x = as.factor(Alpha), y = KLdiv)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape =NA, fill="lightblue") +  
  labs(y = "", x = TeX("$\\alpha$")) +
  theme_minimal() + ggtitle("KL divergence") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) 

ma = ggplot(data, aes(x = as.factor(Alpha), y = MedAPE)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape =NA, fill="lightgreen") +  
  labs(y = "", x = TeX("$\\alpha$")) +
  theme_minimal() + ggtitle("MedAPE") + 
  theme(legend.position = "none" , plot.title = element_text(hjust = 0.5)) +
  ylim(c(0,5))

rr = ggplot(data, aes(x = as.factor(Alpha), y = RRSE)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape =NA, fill="orange") +  
  labs(y = "", x = TeX("$\\alpha$")) +
  theme_minimal() + ggtitle("RRSE") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) + 
  ylim(c(0,4)) 

dev.new(width=12, height=7.5)
(kl | ma | rr)  #boxplots of the indicators depending on the value of alpha


##############  forecasting residuals plots

class = "15_21"
year = "T_24"

to_ggplot_df_singleyear = function(fdata, year, names, label_year = "2024") {
  df = data.frame(cbind(names),fdata$data)
  names(df) = c("name", paste0("val",1:length(fdata$argvals)))
  df = df %>% pivot_longer(cols = paste0("val",1:length(fdata$argvals)), values_to = "dens", names_to = "dummy") %>% dplyr::select(name,dens) 
  df$x = rep(fdata$argvals, 107)
  df$Year = rep(label_year, length(fdata$argvals) * 107)
  df
}

#plot for the alpha-Tsagris residuals
plot_provs_a = to_ggplot_df_singleyear(alpha.folding(residuals[["0.4"]],0.4), year, provinces_names)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
mean_provs = plot_provs_a %>% group_by(x, Year) %>% summarise(mdens = mean(dens))

colors = rep("gray",107)
names(colors) = levels(factor(plot_provs_a$name))
colScale = scale_colour_manual(name = "name", values = colors)
lims = range(as.POSIXct( as.Date(plot_provs_a$x,origin = "1970-01-01")))

tsag = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale +  
  geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + 
  geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = mdens,color = Year), color = "black",size = 1.4)+
  #scale_x_datetime(labels = time_format("%b"),limits = lims) + 
  facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + 
  font("xy.text",size = 17) +
  ylim(c(0.01, 0.036))


#plot for the clr residuals
plot_provs_a = to_ggplot_df_singleyear(inv_clr(residuals[["0"]]), year, provinces_names)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
mean_provs = plot_provs_a %>% group_by(x, Year) %>% summarise(mdens = mean(dens))

lims = range(as.POSIXct( as.Date(plot_provs_a$x,origin = "1970-01-01")))

clr = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale +  
  geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + 
  geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = mdens,color = Year), color = "black",size = 1.4)+
  #scale_x_datetime(labels = time_format("%b"),limits = lims) + 
  facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + 
  font("xy.text",size = 17) + ylim(c(0.01, 0.036))

dev.new(width=12, height=6)
(clr + tsag)


dev.new(width=6, height=12)
(clr / tsag)

