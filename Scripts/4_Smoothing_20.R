rm(list = ls())
graphics.off()
source("Scripts/Utility_Functions.R")
load("Output/Data/Aggregated_Provinces.Rdata") ### for smoothing on Provinces

#### Prepare data format for smoothing
data_to_smooth_provinces = provinces_aggregated %>% tidyr::complete (NOME_PROVINCIA,partial_date_death,CL_ETA)
data_to_smooth_provinces[is.na(data_to_smooth_provinces)] = 0
data_to_smooth_provinces$CL_ETA = factor(data_to_smooth_provinces$CL_ETA, levels = c("0_10","11_14","15_21"))
data_to_smooth_provinces = data_to_smooth_provinces %>% group_by(NOME_PROVINCIA,CL_ETA) %>% dplyr::select(-starts_with("F")) %>%  dplyr::select(-starts_with("M")) %>% arrange(NOME_PROVINCIA)
data_to_smooth_provinces = data_to_smooth_provinces %>% mutate_at(.vars = paste0("T_",11:20),.funs = function(x){return(x/sum(x))})


######### smoothing without 0 imputation

load("Output/Data/Smoothing_datasets.Rdata")

classe = c("0_10","11_14","15_21")
years = paste0("T_",11:20)
alfa = 5e-3
der = 2
ch = 1
npoints = 1000
knots = 52
xcp = as.numeric((as.numeric(as.Date("2020-01-01")):as.numeric(as.Date("2020-12-31"))))
knots_vector = floor(xcp[seq(1,length(xcp), length.out = knots)])
degree = 3
add_node = 1
add_days = 7

xcp_extended = c((xcp[1] - add_days):(xcp[1] - 1),xcp,(xcp[length(xcp)] + 1):(xcp[length(xcp)] + add_days))  #add one week before ad one week after
knots_vector_extended = c(xcp_extended[1],knots_vector,xcp_extended[length(xcp_extended)])  #54 knots
w = rep(1, length(xcp_extended))
start = xcp[1]
stop = xcp[length(xcp)]
argvals = seq(xcp_extended[1],xcp_extended[length(xcp_extended)],length.out = npoints)

alpha.seq = seq(0, 1, by=0.05) #alpha transformation parameter
#the final parameter is to be computed by minimizing the prediction error norm

class = classe[3]  #focus only on elderly class


######## smoothing functions

smoothing.row <- function(alpha, vector) {
  nbasis <- 4:50
  
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
                           knots = knots_vector_extended,num_points = npoints,xcp = xcp_extended, fast = 1)
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


##### TSAGRIS

ncores <- max(parallel::detectCores() -1,1)
cl <- makeCluster(ncores)
clusterExport(cl, 
              c("xcp_extended", "alfa", "w", "ch", "der", "degree", "knots_vector_extended",
                "smoothing.row", "argvals", "xcp_extended","create.fourier.basis", 
                "smooth.basis", "eval.fd", "fdPar","npoints"))
registerDoParallel(cl)


clist_province = flist_province = list()

for (alpha in alpha.seq) {
  
  print("---------------")
  idx = paste0(alpha)
  print(idx)
  clist_province[[idx]] = list()
  
  for(year in years) {
    print(year)
    temp = data_to_smooth_provinces %>% filter(CL_ETA == class) %>% ungroup() %>% dplyr::select(NOME_PROVINCIA,partial_date_death,one_of(year)) %>% arrange(partial_date_death)
    temp$dummy = temp[[year]]
    temp = temp %>% dplyr::select(-one_of(year))
    temp = temp %>% pivot_wider(names_from = partial_date_death, values_from = dummy)
    temp[is.na(temp)] = 0
    matrix = as.matrix(temp[,-1])
    matrix[which(rowSums(matrix) < 1/2),] = 1/dim(matrix)[2]
    matrix = cbind(matrix[,length(xcp):(length(xcp) - add_days + 1)], matrix, matrix[,1:add_days])
    
    clist_province[[idx]][[year]] = smoothing(alpha, matrix)
  }
  
}

for (year in years) flist_province[[year]] = alpha.folding(clist_province[["0.5"]][[year]], 0.5)

save(clist_province,flist_province, file = Rdata_path(paste0("Smoothing_Provinces_",knots,"_alpha")))


#### ISOMETRIC

smoothing <- function(alpha, matrix) {
  
  ret = NULL
  
  if (alpha == 0) {  
    smooth = smoothSplines(k = degree,l = der,alpha = alfa,data = matrix,
                           knots = knots_vector_extended,num_points = npoints,xcp = xcp_extended, fast = 1)
    smoothed = fdata(smooth$Y, argvals = argvals)
    smoothed = conditioning(smoothed, start, stop)
    ret = clr(smoothed)
  }
  
  else {
    
    matrix.alpha = alpha.isometric(fdata(matrix), alpha)$data
    
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


clist_province = flist_province = list()

for (alpha in alpha.seq) {
  
  print("---------------")
  idx = paste0(alpha)
  print(idx)
  clist_province[[idx]] = list()
  
  for(year in years) {
    print(year)
    temp = data_to_smooth_provinces %>% filter(CL_ETA == class) %>% ungroup() %>% dplyr::select(NOME_PROVINCIA,partial_date_death,one_of(year)) %>% arrange(partial_date_death)
    temp$dummy = temp[[year]]
    temp = temp %>% dplyr::select(-one_of(year))
    temp = temp %>% pivot_wider(names_from = partial_date_death, values_from = dummy)
    temp[is.na(temp)] = 0
    matrix = as.matrix(temp[,-1])
    matrix[which(rowSums(matrix) < 1/2),] = 1/dim(matrix)[2]
    matrix = cbind(matrix[,length(xcp):(length(xcp) - add_days + 1)], matrix, matrix[,1:add_days])
    
    clist_province[[idx]][[year]] = smoothing(alpha, matrix)
  }
  
}

for (year in years) flist_province[[year]] = alpha.isometric.inv(clist_province[["0.5"]][[year]], 0.5)

save(clist_province, flist_province, file = Rdata_path(paste0("Smoothing_Provinces_",knots,"_alpha_isometric")))






