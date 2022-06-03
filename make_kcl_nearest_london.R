library(openair)
library(dplyr)
library(worldmet)
library(purrr)
library(tidyr)
library(stringr)

#Read in air quality stations from London source, codes and lat/lon details

london_kcl_meta = importMeta(source = "kcl") %>%
  drop_na()


#read in meteorology stations, codes and lat/lon details

directory_met_data_file = "D:/cpdav/UK_met_data/noaa_UK_met_data_"
noaa_files = list.files("D:/cpdav/UK_met_data/")

noaa_file_names =  noaa_files %>%
  paste0(collapse = " ")

#check data availability
#function data_availability() available in script data_availability_noaa_function.R

noaa_meta_UK = getMeta(country = "UK") %>%
  filter(begin <= as.Date("2016-01-01") & end >= as.Date("2020-12-31"))

noaa_meta_code_list = noaa_meta_UK$code

data_availability_noaa_meta = data_availability(noaa_meta_code_list, directory_met_data_file)

data_availability_noaa_meta_filtered = data_availability_noaa_meta %>%
  filter(wd_na >= 80.)

noaa_meta_filtered_code = data_availability_noaa_meta_filtered$code

noaa_meta_UK = noaa_meta_UK %>%
  filter(code %in% noaa_meta_filtered_code) %>%
  mutate(avaliable = map_lgl(code,~str_detect(noaa_file_names,.x))) %>%
  filter(avaliable)

kcl_sf = sf::st_as_sf(london_kcl_meta,coords = c("longitude","latitude")) %>% 
  mutate(index = 1:nrow(.))

noaa_sf = sf::st_as_sf(noaa_meta_UK,coords = c("longitude","latitude"))

kcl_noaa_nearest = noaa_sf[sf::st_nearest_feature(kcl_sf,noaa_sf),] %>% 
  mutate(index = 1:nrow(.)) %>% 
  group_nest(index, .key = "noaa") %>% 
  left_join(group_nest(kcl_sf,index,.key = "kcl"),"index") %>% 
  mutate(code = map_chr(kcl,~unique(.$code)),
         AQ_site = map_chr(kcl,~unique(.$site)),
         met_code = map_chr(noaa,~unique(.$code)),
         met_site = map_chr(noaa,~unique(.$station))) %>% 
  left_join(london_kcl_meta[,c("code","site_type")],"code") %>%
  filter(met_site %in% c("CITY", "HEATHROW", "KENLEY AIRFIELD", "NORTHOLT"))
