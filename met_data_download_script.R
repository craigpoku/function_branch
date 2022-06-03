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
  filter(begin <= as.Date("2006-01-01") & end >= as.Date("2020-12-31"))

noaa_meta_UK_code = noaa_meta_UK$code

noaa_meta_UK_code_test = noaa_meta_UK_code[0:2]

test = importNOAA(code = "030050-99999", year = 2018:2020)

met_data_download = function(codes, begin_year, end_year) {
  for (i in 1:length(codes)){
    print(c(i, codes[i]))
    met_raw = importNOAA(code = codes[i], year = begin_year:end_year)
    
    write.csv(met_raw, paste("C:/Users/cpdav/Documents/GitHub/UK_met_data_covid/", "noaa_UK_met_data_",toString
                            (codes[i]),".csv"), row.names = FALSE)
  }
}

met_data_download(noaa_meta_UK_code[133:189], 2006, 2020)

