rm(list = ls())
source("Scripts/Utility_Functions.R") #Load packages and functions which are used in the project

########## This chunk of code removes artificial NA (e.g. 29 february for non leap years) and provides better date format #####
modified_date_dataset = read.csv("RawData/comuni_giornaliero_29febbraio24.csv", header = T)
modified_date_dataset$GE = as.character(modified_date_dataset$GE)
months = str_sub(modified_date_dataset$GE,1,-3)
days = str_sub(modified_date_dataset$GE,-2,-1)
modified_date_dataset$partial_date_death = as.Date(paste0("2024-",months,"-",days))
modified_date_dataset$CL_ETA = as.character(modified_date_dataset$CL_ETA)
modified_date_dataset$COD_PROVCOM = as.character(modified_date_dataset$COD_PROVCOM)
modified_date_dataset = modified_date_dataset %>% replace(is.na(.),0) %>% mutate(T_24 = replace(T_24, T_24 == "n.d.", 0), 
                                                                                 M_24 = replace(M_24, M_24 == "n.d.", 0),
                                                                                 F_24 = replace(F_24, F_24 == "n.d.", 0))
modified_date_dataset = modified_date_dataset %>% dplyr::select(-GE)
save(modified_date_dataset, file = Rdata_path("modified_istat_dataset_febbraio24")) 
#### Note: All dates are specified at 2024 for simplicity, but the actual year is specified in each column.
rm(days)
rm(months)


#### Aggragation by age (See notation in the pdf from ISTAT) ######
firstvec = c(0,11,15)
lastvec = c(10,14,21)
#### So we are aggregating 0-49 years, 50-69 years, 70+ years #####

provinces_names = sort(unique(modified_date_dataset$NOME_PROVINCIA))
modified_date_dataset = modified_date_dataset %>% mutate_at(c('T_24', 'M_24','F_24'), as.numeric)
provinces_data = modified_date_dataset %>% 
  group_by(NOME_PROVINCIA,CL_ETA,partial_date_death) %>% 
  summarise_at(c(paste0("M_",c(11:24)),paste0("F_",c(11:24)),paste0("T_",c(11:24))), sum) 
provinces_aggregated = aggregate_classes(provinces_data,keyname = "NOME_PROVINCIA",firstvec,lastvec)
save(list = c("provinces_aggregated","provinces_names"), file = Rdata_path("Aggregated_Provinces_24"))



#### We have to adjust same names of the provinces
load("Output/Data/Aggregated_Provinces_24.Rdata")
source("Scripts/Utility_Functions.R")
provinces_data = read.csv("RawData/Provinces_Population_2023.csv")
provinces_data = provinces_data %>% dplyr::select(Territorio, Sesso, Età, Value) 
provinces_data$Territorio[provinces_data$Territorio =="Valle d'Aosta / Vallée d'Aoste"] = "Valle d'Aosta/Vall\xe9e d'Aoste"
provinces_data$Territorio[provinces_data$Territorio =="Bolzano / Bozen"] ="Bolzano/Bozen"
provinces_data$Territorio[provinces_data$Territorio =="Reggio di Calabria"] = "Reggio Calabria"
provinces_data$Territorio[provinces_data$Territorio =="Forlì-Cesena"] = "Forl\xec-Cesena"
provinces_data = provinces_data %>% filter(Territorio %in% provinces_names, Sesso =="totale")
provinces_data = distinct(provinces_data)
provinces_data$Età =  substr(provinces_data$Età,start=1, stop= nchar(provinces_data$Età) - 5)
provinces_data$Età[provinces_data$Età == "100 anni "] = "100"
provinces_data$Età = as.numeric(provinces_data$Età) 
provinces_data = na.omit(provinces_data)
provinces_data = provinces_data %>% mutate(CL_ETA = age_to_fact(Età)) %>% 
  dplyr::select(-c(Età,Sesso)) %>% group_by(Territorio,CL_ETA) %>% summarise(num_residenti = sum(Value))
colnames(provinces_data) = c("Provincia", "CL_ETA", "num_residenti")

save(provinces_aggregated,provinces_names, provinces_data, file = Rdata_path("Aggregated_Provinces_24"))


