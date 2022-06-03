data_availability = function(codes, directory_met_data_file) {
  code = as.character(codes)
  site_name = as.character(codes)
  wd_na = as.numeric(codes)
  ws_na = as.numeric(codes)
  air_temp_na = as.numeric(codes)
  atmos_pres_na = as.numeric(codes)

  
  for (i in 1:length(codes)){
    x = paste(directory_met_data_file,toString
              (codes[i]),".csv", sep = "")
    print(c(i, codes[i], x, file.size(x)))
    if (file.exists(x) & file.size(x) >= 50) {
      file_i = read.csv(x)
 
      code[i] = codes[i]
      site_name[i] = unique(file_i$station)
      wd_na[i] = 100-(sum(is.na(file_i$wd))/nrow(file_i))*100
      ws_na[i] = 100-(sum(is.na(file_i$ws))/nrow(file_i))*100
      air_temp_na[i] = 100-(sum(is.na(file_i$air_temp))/nrow(file_i))*100
      atmos_pres_na[i] = 100-(sum(is.na(file_i$atmos_pres))/nrow(file_i))*100
      
    } 
  }
  
  data.frame(code, site_name, wd_na, ws_na, air_temp_na, atmos_pres_na, stringsAsFactors=FALSE)
}
  
