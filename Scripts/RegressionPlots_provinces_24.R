######################################################## Regression plots ##############
rm(list = ls())
source("Scripts/Utility_Functions.R")
load("Output/Data/Aggregated_Provinces_24.Rdata")
load("Output/Data/Predictive_performance_24.Rdata")

classe = c("0_10","11_14","15_21")
years = paste0("T_",11:24)
class_labs = c("0-49 years", "50-69 years", "70+ years")
year_labs = paste0("20",11:24)
names(class_labs) =  classe
names(year_labs) = years
folder_names = c("0-49","50-69","70+")
names(folder_names) = classe
years_to_plot = years[14]

class = "15_21"
year = "T_24"

alpha.best = 0.85
alpha.best_isom = 1

temp_res_clr = inv_clr(clr_performance$errors)
temp_res = alpha.folding(tsagris_performance$errors, alpha.best)
temp_res_isom = alpha.isometric.inv(isometric_performance$errors, alpha.best_isom)

to_ggplot_df_singleyear = function(fdata, year, names, label_year = "2024") {
  df = data.frame(cbind(names),fdata$data)
  names(df) = c("name", paste0("val",1:length(fdata$argvals)))
  df = df %>% pivot_longer(cols = paste0("val",1:length(fdata$argvals)), values_to = "dens", names_to = "dummy") %>% dplyr::select(name,dens) 
  df$x = rep(fdata$argvals, 107)
  df$Year = rep(label_year, length(fdata$argvals) * 107)
  df
}

#plot for the alpha-Tsagris residuals
plot_provs_a = to_ggplot_df_singleyear(temp_res, year, provinces_names)

colors = rep("gray",107)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))

plot_provs_alpha = to_ggplot_df_singleyear(tsagris_performance$errors,year, provinces_names )
plot_provs_alpha$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
provinces_data_t = provinces_data %>% filter(CL_ETA == class)
provinces_data_t = provinces_data_t %>% mutate(weight = num_residenti/sum(provinces_data$num_residenti)) %>% dplyr::select(Provincia, weight)
plot_provs_alpha = plot_provs_alpha %>% left_join(provinces_data_t, by = c("name" = "Provincia"))
plot_provs_alpha =  plot_provs_alpha %>% mutate(cl_dens = dens*weight)
mean_provs = plot_provs_alpha %>% group_by(x, Year) %>% summarise(mdens = mean(dens))
Integrals = mean_provs %>% group_by(Year) %>% summarise(Integrals = int.simpson2(x = x,y = (mdens*alpha.best + 1)^(1/alpha.best)))
mean_provs = mean_provs %>% left_join(Integrals)
mean_provs = mean_provs %>% group_by(Year,x) %>% mutate(meandens = 
                                                          (mdens*alpha.best + 1)^(1/alpha.best)/Integrals)
names(colors) = levels(factor(plot_provs_a$name))
colScale = scale_colour_manual(name = "name", values = colors)
lims = range(as.POSIXct( as.Date(plot_provs_a$x,origin = "1970-01-01")))

variance = plot_provs_a %>% group_by(name) %>% summarise(v = var(dens)) 
v = mean(variance$v)

g = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale +   geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = meandens,color = Year), color = "black",size = 1.4)+
  #scale_x_datetime(labels = time_format("%b"),limits = lims) + 
  facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + font("xy.text",size = 17) #+ ggtitle(paste("Mean variance of the curves:",v))
#x11()
pdf(paste0("Output/Plot/",folder_names[[class]],"/Res_regr_Tsagris_24.pdf"), width = 20, height = 10)

plot(g)
dev.off()



#plot for the clr residuals
plot_provs_a = to_ggplot_df_singleyear(temp_res_clr, year, provinces_names)

colors = rep("gray",107)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))

plot_provs_alpha = to_ggplot_df_singleyear(clr_performance$errors, year, provinces_names )
plot_provs_alpha$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
provinces_data_t = provinces_data %>% filter(CL_ETA == class)
provinces_data_t = provinces_data_t %>% mutate(weight = num_residenti/sum(provinces_data$num_residenti)) %>% dplyr::select(Provincia, weight)
plot_provs_alpha = plot_provs_alpha %>% left_join(provinces_data_t, by = c("name" = "Provincia"))
plot_provs_alpha =  plot_provs_alpha %>% mutate(cl_dens = dens*weight)
mean_provs = plot_provs_alpha %>% group_by(x, Year) %>% summarise(mdens = mean(dens))
Integrals = mean_provs %>% group_by(Year) %>% summarise(Integrals = int.simpson2(x = x,y = exp(mdens)))
mean_provs = mean_provs %>% left_join(Integrals)
mean_provs = mean_provs %>% group_by(Year,x) %>% mutate(meandens = exp(mdens)/Integrals)
names(colors) = levels(factor(plot_provs_a$name))
colScale = scale_colour_manual(name = "name", values = colors)
lims = range(as.POSIXct( as.Date(plot_provs_a$x,origin = "1970-01-01")))

variance = plot_provs_a %>% group_by(name) %>% summarise(v = var(dens)) 
v = mean(variance$v)

g = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale +   geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = meandens,color = Year), color = "black",size = 1.4)+
  #scale_x_datetime(labels = time_format("%b"),limits = lims) + 
  facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + font("xy.text",size = 17) #+ ggtitle(paste("Mean variance of the curves:",v))
#x11()
pdf(paste0("Output/Plot/",folder_names[[class]],"/Res_regr_clr_24.pdf"), width = 20, height = 10)

plot(g)
dev.off()


## plot for the isometric residuals

plot_provs_a = to_ggplot_df_singleyear(temp_res_isom, year, provinces_names)

colors = rep("gray",107)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))

plot_provs_alpha = to_ggplot_df_singleyear(isometric_performance$errors,year, provinces_names )
plot_provs_alpha$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
provinces_data_t = provinces_data %>% filter(CL_ETA == class)
provinces_data_t = provinces_data_t %>% mutate(weight = num_residenti/sum(provinces_data$num_residenti)) %>% dplyr::select(Provincia, weight)
plot_provs_alpha = plot_provs_alpha %>% left_join(provinces_data_t, by = c("name" = "Provincia"))
plot_provs_alpha =  plot_provs_alpha %>% mutate(cl_dens = dens*weight)
mean_provs = plot_provs_alpha %>% group_by(x, Year) %>% summarise(mdens = mean(dens))

mean_provs$meandens = mean_provs$mdens #initialization

a = matrix(nrow=length(unique(mean_provs$Year)), ncol=length(unique(mean_provs$x)))
for (i in 1:length(unique(mean_provs$Year))) {
  yy = unique(mean_provs$Year)[i]
  v = mean_provs %>% filter(Year == yy) %>% dplyr::select(mdens)
  a[i,] = alpha.isometric.inv.vec(x=v$x, y=v$mdens, alpha=alpha.best_isom)
}

rownames(a) = unique(mean_provs$Year)
colnames(a) = unique(mean_provs$x)

for (i in 1:length(unique(mean_provs$Year))) {
  
  yy = unique(mean_provs$Year)[i]
  
  for (j in 1:length(unique(mean_provs$x))) {
    
    xx = unique(mean_provs$x)[j]
    mean_provs[ mean_provs$x == xx & mean_provs$Year == yy, 4] = a[i,j]
    
  }
}


names(colors) = levels(factor(plot_provs_a$name))
colScale = scale_colour_manual(name = "name", values = colors)
lims = range(as.POSIXct( as.Date(plot_provs_a$x,origin = "1970-01-01")))

variance = plot_provs_a %>% group_by(name) %>% summarise(v = var(dens)) 
v = mean(variance$v)

g = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale +   geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = meandens,color = Year), color = "black",size = 1.4)+
  #scale_x_datetime(labels = time_format("%b"),limits = lims) + 
  facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + font("xy.text",size = 17) #+ ggtitle(paste("Mean variance of the curves:",v))
#x11()
pdf(paste0("Output/Plot/",folder_names[[class]],"/Res_regr_Isometric_24.pdf"), width = 20, height = 10)

plot(g)
dev.off()



