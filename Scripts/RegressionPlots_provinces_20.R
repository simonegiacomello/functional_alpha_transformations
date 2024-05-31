######################################################## Regression plots ##############
rm(list = ls())
source("Scripts/Utility_Functions.R")
load("Output/Data/Aggregated_Provinces.Rdata")
load("Output/Data/Regression_Province_clr.Rdata")
load("Output/Data/Regression_Province_best_alpha.Rdata")
load("Output/Data/Regression_Province_best_alpha_isometric.Rdata")


classe = c("0_10","11_14","15_21")
years = paste0("T_",11:20)
class_labs = c("0-49 years", "50-69 years", "70+ years")
year_labs = paste0("20",11:20)
names(class_labs) =  classe
names(year_labs) = years
folder_names = c("0-49","50-69","70+")
names(folder_names) = classe
years_to_plot = years[7:10]
class = classe[3]

temp_res_clr = errors_province_clr
temp_res = errors_province_best
temp_res_isom = errors_province_best_isom

for(year in years_to_plot) {
  temp_res_clr[[year]] = inv_clr(temp_res_clr[[year]])
  temp_res[[year]] = alpha.folding(temp_res[[year]], alpha.best)
  temp_res_isom[[year]] = alpha.isometric.inv(temp_res_isom[[year]], alpha.best_isom)
}


#plot for the alpha-Tsagris residuals

plot_provs_a = fdata_to_ggplot_df(fdatalist = temp_res,names = provinces_names, years = c("T_17","T_18","T_19","T_20"), year_labs = year_labs[7:10] )
colors = rep("gray",107)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
plot_provs_alpha = fdata_to_ggplot_df(fdatalist = errors_province_best,names = provinces_names, years = c("T_17","T_18","T_19","T_20"), year_labs = year_labs[7:10] )
plot_provs_alpha$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
provinces_data_t = provinces_data %>% filter(CL_ETA == class)
provinces_data_t = provinces_data_t %>% mutate(weight = num_residenti/sum(provinces_data$num_residenti)) %>% dplyr::select(provincia, weight)
plot_provs_alpha = plot_provs_alpha %>% left_join(provinces_data_t, by = c("name" = "provincia"))
plot_provs_alpha =  plot_provs_alpha %>% mutate(cl_dens = dens*weight)
mean_provs = plot_provs_alpha %>% group_by(x, Year) %>% summarise(mdens = mean(dens))
Integrals = mean_provs %>% group_by(Year) %>% summarise(Integrals = int.simpson2(x = x,y = (mdens*alpha.best + 1)^(1/alpha.best)))
mean_provs = mean_provs %>% left_join(Integrals)
mean_provs = mean_provs %>% group_by(Year,x) %>% mutate(meandens =
                                                          (mdens*alpha.best + 1)^(1/alpha.best)/Integrals)
names(colors) = levels(factor(plot_provs_a$name))
colScale = scale_colour_manual(name = "name", values = colors)
lims = range(as.POSIXct( as.Date(plot_provs_a$x,origin = "1970-01-01")))

# v = var(mean_provs$meandens)
variance = plot_provs_a %>% group_by(name) %>% summarise(v = var(dens))
v = mean(variance$v)

g = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale  +  ylim(0,0.02) +   geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = meandens,color = Year), color = "black",size = 1.4)+
  scale_x_datetime(labels = time_format("%b"),limits = lims) + facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + font("xy.text",size = 17) + ggtitle(paste("Mean variance of the curves:",v))
#x11()
pdf(paste0("Output/Plot/",folder_names[[class]],"/Res_regr_Tsagris_20.pdf"), width = 20, height = 10)

plot(g)
dev.off()



#plot for the clr residuals
plot_provs_c = fdata_to_ggplot_df(fdatalist = temp_res_clr,names = provinces_names, years = c("T_17","T_18","T_19","T_20"), year_labs = year_labs[7:10] )
colors = rep("gray",107)
plot_provs_c$name = factor(plot_provs_c$name, levels = c(setdiff(plot_provs_c$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
plot_provs_clr = fdata_to_ggplot_df(fdatalist = errors_province_clr,names = provinces_names, years = c("T_17","T_18","T_19","T_20"), year_labs = year_labs[7:10] )
plot_provs_clr$name = factor(plot_provs_c$name, levels = c(setdiff(plot_provs_c$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
provinces_data_t = provinces_data %>% filter(CL_ETA == class)
provinces_data_t = provinces_data_t %>% mutate(weight = num_residenti/sum(provinces_data$num_residenti)) %>% dplyr::select(provincia, weight)
plot_provs_clr = plot_provs_clr %>% left_join(provinces_data_t, by = c("name" = "provincia"))
plot_provs_clr =  plot_provs_clr %>% mutate(cl_dens = dens*weight)
mean_provs = plot_provs_clr %>% group_by(x, Year) %>% summarise(mdens = mean(dens))
Integrals = mean_provs %>% group_by(Year) %>% summarise(Integrals = int.simpson2(x = x,y = exp(mdens)))
mean_provs = mean_provs %>% left_join(Integrals)
mean_provs = mean_provs %>% group_by(Year,x) %>% mutate(meandens = exp(mdens)/Integrals)
names(colors) = levels(factor(plot_provs_c$name))
colScale = scale_colour_manual(name = "name", values = colors)
lims = range(as.POSIXct( as.Date(plot_provs_c$x,origin = "1970-01-01")))

#v = var(mean_provs$meandens)
variance = plot_provs_c %>% group_by(name) %>% summarise(v = var(dens))
v = mean(variance$v)

g = ggplot(plot_provs_c, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale  +  ylim(0,0.02) +   geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  + geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = meandens,color = Year), color = "black",size = 1.4)+
  scale_x_datetime(labels = time_format("%b"),limits = lims) + facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + font("xy.text",size = 17) + ggtitle(paste("Mean variance of the curves:",v))
#x11()
pdf(paste0("Output/Plot/",folder_names[[class]],"/Res_regr_clr_20.pdf"), width = 20, height = 10)

plot(g)
dev.off()


#plot for the alpha-isometric residuals

plot_provs_a = fdata_to_ggplot_df(fdatalist = temp_res_isom,names = provinces_names, years = c("T_17","T_18","T_19","T_20"), year_labs = year_labs[7:10] )
colors = rep("gray",107)
plot_provs_a$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
plot_provs_alpha = fdata_to_ggplot_df(fdatalist = errors_province_best_isom,names = provinces_names, years = c("T_17","T_18","T_19","T_20"), year_labs = year_labs[7:10] )
plot_provs_alpha$name = factor(plot_provs_a$name, levels = c(setdiff(plot_provs_a$name, c("Roma", "Milano", "Bergamo","Napoli") ), c("Roma", "Milano", "Bergamo","Napoli")))
provinces_data_t = provinces_data %>% filter(CL_ETA == class)
provinces_data_t = provinces_data_t %>% mutate(weight = num_residenti/sum(provinces_data$num_residenti)) %>% dplyr::select(provincia, weight)
plot_provs_alpha = plot_provs_alpha %>% left_join(provinces_data_t, by = c("name" = "provincia"))
plot_provs_alpha =  plot_provs_alpha %>% mutate(cl_dens = dens*weight)
mean_provs = plot_provs_alpha %>% group_by(x, Year) %>% summarise(mdens = mean(dens))
mean_provs$meandens = mean_provs$mdens  #initialization

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

# v = var(mean_provs$meandens)
variance = plot_provs_a %>% group_by(name) %>% summarise(v = var(dens)) 
v = mean(variance$v) 

g = ggplot(plot_provs_a, mapping = aes(x = as.POSIXct( as.Date(x,origin = "1970-01-01")), y = dens, color = name) )+ colScale  +  ylim(0,0.02) +   geom_line(aes(alpha  = ifelse(name %in% c("Roma", "Milano", "Bergamo","Napoli"), 1, 1)), size = 1.1 )  +
  geom_line(data = mean_provs, mapping=aes(x = as.POSIXct(as.Date(x,origin = "1970-01-01")), y = meandens, color = Year), color = "black",size = 1.4)+
  scale_x_datetime(labels = time_format("%b"),limits = lims) + facet_grid(~Year) + theme_pubr(base_size = 20)  + rremove("legend")+ rremove("xlab") + rremove("ylab")  + labs_pubr(base_size = 20) + font("xy.text",size = 17) + ggtitle(paste("Mean variance of the curves:",v))
#x11()
pdf(paste0("Output/Plot/",folder_names[[class]],"/Res_regr_Isometric_20.pdf"), width = 20, height = 10)

plot(g)
dev.off()








