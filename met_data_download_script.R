library(dplyr)
library(rmweather)
library(ranger)
library(ggplot2)
library(worldmet)
library(openair)
library(tidyr)
library(stringr)
library(MAP)
library(purrr)


noaa_files = list.files("C:/Users/cpdav/Documents/GitHub/UK_met_data/")

noaa_file_names =  noaa_files %>%
  paste0(collapse = " ")

noaa_meta_UK = getMeta(country = "UK") %>%
  filter(begin <= as.Date("2006-01-01") & end >= as.Date("2021-12-31"))

noaa_meta_UK_code = noaa_meta_UK$code

test = importNOAA(code = "030050-99999", year = 2018:2020)

met_data_download = function(codes, begin_year, end_year) {
  for (i in 1:length(codes)){
    print(c(i, codes[i]))
    met_raw = importNOAA(code = codes[i], year = begin_year:end_year)
    
    write.csv(met_raw, paste("/home/craig/UK_met_data/", "noaa_UK_met_data_",toString
                            (codes[i]),".csv"), row.names = FALSE, sep = "")
  }
}

met_data_download(noaa_meta_UK_code, 2011, 2021)

