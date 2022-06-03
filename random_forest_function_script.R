#Functions specific to running random forest code for sites across UK

#----------- Creates dataset that combines AQ and met dataframes, whilst removing any non-existing data, 5 variables --------------------


read_met_sites = function(site_code, metcode, directory_met_data, begin_date, end_date){
  
  met_UK_df = read.csv(paste(directory_met_data, metcode,".csv")) %>%
    mutate(date = lubridate::ymd_hms(date)) %>%
    filter(date >= as.Date(begin_date) & date <= as.Date(end_date)) %>%
    rename(met_code = code)
  
  begin_year = as.numeric(substring(begin_date,1,4))
  end_year = as.numeric(substring(end_date,1,4))
  
  UK_raw_df = importAURN(site = site_code, year = begin_year:end_year, meta = TRUE)
  
  if(!is.null(UK_raw_df)){
    met_UK_df = UK_raw_df %>% 
      select(-any_of(c("ws", "wd", "air_temp", "latitude", "longitude"))) %>% 
      left_join(met_UK_df, ., by = "date")
    print(paste0(site_code, "/", metcode, " finished."))
  } else {
    print(paste0(site_code, "/", metcode, " is missing."))
  }
  
  return(tibble::tibble(met_UK_df))
  
}

#--------------- Using combined AQ_met dataframe, function prepares dataset to be suitable for rmweather  ----------------------

met_aq_prepared_rm = function(df, pollutant_type, site_type_category){
  
  df_rm = df %>%
    select(code, site_type, met_code, station, date, all_of(pollutant_type), 
           ws, wd, air_temp, atmos_pres, RH, cl)%>%
    filter(site_type == site_type_category) %>%
    pivot_longer(-c(code, site_type, met_code, station, date, ws, wd, air_temp, atmos_pres, RH, cl),
                 names_to = "pollutant") %>%
    select(-site_type) %>%
    rename(atmospheric_pressure = atmos_pres) %>%
    filter(!is.na(value))
  
}

#--------------Creates df that contains observed vs BAU data, as well as normalised concentration  for all sites specified in the prepared rmweather df--------------
#also contains relevant statistics e.g. MSE, rsquared coefficient used to check viability and bias in output

normalised_concentration_and_BAU_observed_combined_df = function(df, site, n_tree, begin_date_train, end_date_train, 
                                                              end_date_predict, data_threshold, n_sample){
  
  #1. checks there is enough data for random forest to work
  
  print(paste(site, "rm process beginning"))
  
  begin_year = as.numeric(substring(begin_date_train,1,4))
  end_year = as.numeric(substring(end_date_train,1,4))
  
  begin_month = as.numeric(substring(begin_date_train,6,7))
  end_month = as.numeric(substring(end_date_train,6,7))
  
  begin_day = as.numeric(substring(begin_date_train,9,10))
  end_day = as.numeric(substring(end_date_train,9,10))  
  
  diff_hours = c(ISOdate(end_year,end_month,end_day), ISOdate(begin_year,begin_month,begin_day))
  
  diff_hours_quality_check = difftime(diff_hours[1], diff_hours[2], units="hours")+1
  
  df_quality_check = df %>% 
    filter(date >= as.Date(begin_date_train), code == site)
  
  total_date_df = length(df_quality_check$date)
  print(length(df_quality_check$date)/as.numeric(diff_hours_quality_check))
  
  if(length(df_quality_check$date)/as.numeric(diff_hours_quality_check) >= data_threshold){
    
    df_prepared_train = df_quality_check %>% 
      filter(date <= as.Date(end_date_train))%>%
      rmw_prepare_data(na.rm = TRUE)
    
    
    df_prepared_predict = df_quality_check %>% 
      filter(date <= as.Date(end_date_predict)) %>%
      rmw_prepare_data(na.rm = TRUE)
    
    
    variables_atmos = c("ws","wd","air_temp",
                        "RH", "atmospheric_pressure", "date_unix", "day_julian", "weekday", "hour")
    
    variables_no_atmos = c("ws","wd","air_temp",
                           "RH", "date_unix", "day_julian", "weekday", "hour")
    
    
    if(sum(df_prepared_predict$atmospheric_pressure, na.rm=T)>0){
      
      #2. performs random forest algorithm using rmweather, both for normalisation and BAU
      
      print(paste(site, "df contains atmospheric pressure variable"))
      rm_df_train = rmw_train_model(df_prepared_train,
                                    variables = variables_atmos,
                                    n_trees = n_tree,
                                    verbose = TRUE)
      
      rm_df_normal = rmw_do_all(df_prepared_predict,
                                variables = variables_atmos,
                                n_trees = n_tree,
                                n_samples = n_tree,
                                verbose = TRUE)
      
      rm_df_train_test = rmw_predict_the_test_set(
        rm_df_train, 
        df = df_prepared_train
      )
      
      testing_mse = as.numeric(sum((rm_df_train_test$value-rm_df_train_test$value_predict)^2)/
                                 count(rm_df_train_test))
      testing_rsquared = cor(rm_df_train_test$value, 
                             rm_df_train_test$value_predict)^2
      
      
      #3. if we wanted to change BAU at this point, we could
      
      rm_df_predict = df_prepared_predict %>% mutate(value_predict = rmw_predict(
        rm_df_train, 
        df = df_prepared_predict),
        rm_normal_value = rm_df_normal$normalised$value_predict,
        training_rsquared = rm_df_train$r.squared, 
        training_mse = rm_df_train$prediction.error,
        testing_rsquared = testing_rsquared, testing_mse = testing_mse) %>%
        select(date, value, rm_normal_value, value_predict, training_rsquared, training_mse, testing_rsquared,
               testing_mse)
      
    } else {
      
      print(paste(site, "df does not contains atmospheric pressure variable"))
      rm_df_train = rmw_train_model(df_prepared_train,
                                    variables = variables_no_atmos,
                                    n_trees = n_tree,
                                    verbose = TRUE)
      
      rm_df_normal = rmw_do_all(df_prepared_predict,
                                variables = variables_no_atmos,
                                n_trees = n_tree,
                                n_samples = n_tree,
                                verbose = TRUE)
      
      rm_df_train_test = rmw_predict_the_test_set(
        rm_df_train, 
        df = df_prepared_train
      )
      
      testing_mse = as.numeric(sum((rm_df_train_test$value-rm_df_train_test$value_predict)^2)/
                                 count(rm_df_train_test))
      testing_rsquared = cor(rm_df_train_test$value, 
                             rm_df_train_test$value_predict)^2
      
      
      rm_df_predict = df_prepared_predict %>% mutate(value_predict = rmw_predict(
        rm_df_train, 
        df = df_prepared_predict),
        rm_normal_value = rm_df_normal$normalised$value_predict,
        training_rsquared = rm_df_train$r.squared, 
        training_mse = rm_df_train$prediction.error,
        testing_rsquared = testing_rsquared, testing_mse = testing_mse) %>%
        select(date, value, rm_normal_value, value_predict, training_rsquared, training_mse, testing_rsquared,
               testing_mse)
      
      
      
      print(paste(site, "rm process complete"))
      print("---------------------------------------------------------------------")
      return(rm_df_predict)
    }}
  else{print(paste(site, "skipped due to insufficient data"))
    print("---------------------------------------------------------------------")
  }
}

#---------------------------Creates df that contains observed vs BAU data using rmweather for all sites specified in the prepared rmweather df --------------
#also contains relevant statistics e.g. MSE, rsquared coefficient for linear regression
# note, there is no normalisation data in this function

BAU_observed_combined_df = function(df, site, n_tree, begin_date_train, end_date_train, 
                                    end_date_predict, data_threshold){
  
  #1. checks there is enough data for random forest to work
  
  print(paste(site, "rm process beginning"))
  
  begin_year = as.numeric(substring(begin_date_train,1,4))
  end_year = as.numeric(substring(end_date_train,1,4))
  
  begin_month = as.numeric(substring(begin_date_train,6,7))
  end_month = as.numeric(substring(end_date_train,6,7))
  
  begin_day = as.numeric(substring(begin_date_train,9,10))
  end_day = as.numeric(substring(end_date_train,9,10))  
  
  diff_hours = c(ISOdate(end_year,end_month,end_day), ISOdate(begin_year,begin_month,begin_day))
  
  diff_hours_quality_check = difftime(diff_hours[1], diff_hours[2], units="hours")+1
  
  df_quality_check = df %>% 
    filter(date >= as.Date(begin_date_train), code == site)
  
  total_date_df = length(df_quality_check$date)
  print(length(df_quality_check$date)/as.numeric(diff_hours_quality_check))
  
  if(length(df_quality_check$date)/as.numeric(diff_hours_quality_check) >= data_threshold){
    
    #2. performs random forest algorithm using rmweather, both for normalisation and BAU
    
    df_prepared_train = df_quality_check %>% 
      filter(date <= as.Date(end_date_train))%>%
      rmw_prepare_data(na.rm = TRUE)
    
    
    df_prepared_predict = df_quality_check %>% 
      filter(date <= as.Date(end_date_predict)) %>%
      rmw_prepare_data(na.rm = TRUE)
    
    
    variables_atmos = c("ws","wd","air_temp",
                        "RH", "atmospheric_pressure", "date_unix", "day_julian", "weekday", "hour")
    
    variables_no_atmos = c("ws","wd","air_temp",
                           "RH", "date_unix", "day_julian", "weekday", "hour")
    
    
    if(sum(df_prepared_predict$atmospheric_pressure, na.rm=T)>0){
      
      print(paste(site, "df contains atmospheric pressure variable"))
      rm_df_train = rmw_train_model(df_prepared_train,
                                    variables = variables_atmos,
                                    n_trees = n_tree,
                                    verbose = TRUE)
      
      
      rm_df_train_test = rmw_predict_the_test_set(
        rm_df_train, 
        df = df_prepared_train
      )
      
      testing_mse = as.numeric(sum((rm_df_train_test$value-rm_df_train_test$value_predict)^2)/
                                 count(rm_df_train_test))
      testing_rsquared = cor(rm_df_train_test$value, 
                             rm_df_train_test$value_predict)^2
      
      #3. if we wanted to change BAU at this point, we could
      
      rm_df_predict = df_prepared_predict %>% mutate(value_predict = rmw_predict(
        rm_df_train, 
        df = df_prepared_predict),
        training_rsquared = rm_df_train$r.squared, 
        training_mse = rm_df_train$prediction.error,
        testing_rsquared = testing_rsquared, testing_mse = testing_mse) %>%
        select(date, value, value_predict, training_rsquared, training_mse, testing_rsquared,
               testing_mse)
      
    } else {
      
      print(paste(site, "df does not contains atmospheric pressure variable"))
      rm_df_train = rmw_train_model(df_prepared_train,
                                    variables = variables_no_atmos,
                                    n_trees = n_tree,
                                    verbose = TRUE)
      
      rm_df_train_test = rmw_predict_the_test_set(
        rm_df_train, 
        df = df_prepared_train
      )
      
      testing_mse = as.numeric(sum((rm_df_train_test$value-rm_df_train_test$value_predict)^2)/
                                 count(rm_df_train_test))
      testing_rsquared = cor(rm_df_train_test$value, 
                             rm_df_train_test$value_predict)^2
      
      
      rm_df_predict = df_prepared_predict %>% mutate(value_predict = rmw_predict(
        rm_df_train, 
        df = df_prepared_predict),
        training_rsquared = rm_df_train$r.squared, 
        training_mse = rm_df_train$prediction.error,
        testing_rsquared = testing_rsquared, testing_mse = testing_mse) %>%
        select(date, value, value_predict, training_rsquared, training_mse, testing_rsquared,
               testing_mse)
      
      
      
      print(paste(site, "rm process complete"))
      print("---------------------------------------------------------------------")
      return(rm_df_predict)
    }}
  else{print(paste(site, "skipped due to insufficient data"))
    print("---------------------------------------------------------------------")
  }
}

#---------------reformats output from random forest df to make more user friendly, add in statistics and has the option to include just a subset of sites (e.g. ULEZ vs Greater London) ------
reformat_random_forest_df_output_statistics_addeded = function(df, UK_code, site_subset_code,
                                                 normal=TRUE, site_subset = TRUE){
  names(df) = UK_code
  
  if(normal==TRUE){
    if(site_subset==TRUE){
      predict_df = plyr::ldply(df, 
                               data.frame) %>%
        filter(!is.na(date)) %>%
        rename(sites = .id, normal_value = rm_normal_value, observed = raw_value)%>%
        select(sites, date, normal_value, observed)%>%
        filter(sites %in% site_subset_code) 
    }
    else{
      predict_df = plyr::ldply(df, 
                               data.frame) %>%
        filter(!is.na(date)) %>%
        rename(sites = .id, normal_value = rm_normal_value, observed = raw_value)%>%
        select(sites, date, normal_value, observed) 
    }
    
    
    predict_df_mean = predict_df %>%
      timeAverage(avg.time = "1 day", statistic = c("mean"))%>%
      rename(normal_mean = normal_value, observed_mean = observed)
    
    predict_df_sd = predict_df %>%
      timeAverage(avg.time = "1 day", statistic = c("sd")) %>%
      rename(normal_sd = normal_value, observed_sd = observed)
    
    predict_df_frequency = predict_df %>%
      select(-normal_value) %>%
      timeAverage(avg.time = "1 day", statistic = c("frequency"))%>%
      rename(observed_frequency = observed)
    
    predict_df_combo = left_join(predict_df_mean, predict_df_sd, by = "date") %>% 
      left_join(., predict_df_frequency, by = "date")
    
    predict_df_combo = predict_df_combo %>%
      mutate(CI_lower = observed_mean - 1.96*(observed_sd/sqrt(observed_frequency)),
             CI_upper = observed_mean + 1.96*(observed_sd/sqrt(observed_frequency))) %>%
      select(-observed_frequency) %>%
      mutate(d7_rollavg_normal_mean = zoo::rollapply(normal_mean,7,mean,align='right',fill=NA),
             d7_rollavg_normal_sd = zoo::rollapply(normal_sd,7,mean,align='right',fill=NA),
             d7_rollavg_observed_mean = zoo::rollapply(observed_mean,7,mean,align='right',fill=NA),
             d7_rollavg_observed_sd = zoo::rollapply(observed_sd,7,mean,align='right',fill=NA),
             d7_rollavg_CI_lower = zoo::rollapply(CI_lower,7,mean,align='right',fill=NA),
             d7_rollavg_CI_upper = zoo::rollapply(CI_upper,7,mean,align='right',fill=NA))
    
    return(predict_df_combo)
  }
  else{
    if(site_subset==TRUE){
      predict_df = plyr::ldply(df, 
                               data.frame) %>%
        filter(!is.na(date)) %>%
        rename(sites = .id, BAU = value_predict, observed = value)%>%
        select(sites, date, BAU, observed)%>%
        filter(sites %in% site_subset_code) 
    }
    else{
      predict_df = plyr::ldply(df, 
                               data.frame) %>%
        filter(!is.na(date)) %>%
        rename(sites = .id, BAU = value_predict, observed = value)%>%
        select(sites, date, BAU, observed)  
    }
    predict_df_mean = predict_df %>%
      timeAverage(avg.time = "1 day", statistic = c("mean"))%>%
      rename(BAU_mean = BAU, observed_mean = observed) %>%
      mutate(delta_BAU_predict_mean = observed_mean- BAU_mean)
    
    predict_df_sd = predict_df %>%
      timeAverage(avg.time = "1 day", statistic = c("sd")) %>%
      rename(BAU_sd = BAU, observed_sd = observed) %>%
      mutate(delta_BAU_predict_sd = abs(observed_sd-BAU_sd))
    
    predict_df_frequency = predict_df %>%
      select(-BAU) %>%
      timeAverage(avg.time = "1 day", statistic = c("frequency"))%>%
      rename(observed_frequency = observed)
    
    predict_df_combo = left_join(predict_df_mean, predict_df_sd, by = "date") %>% 
      left_join(., predict_df_frequency, by = "date")
    
    predict_df_combo = predict_df_combo %>%
      mutate(CI_lower = observed_mean - 1.96*(observed_sd/sqrt(observed_frequency)),
             CI_upper = observed_mean + 1.96*(observed_sd/sqrt(observed_frequency))) %>%
      select(-observed_frequency) %>%
      mutate(d7_rollavg_delta_BAU_predict_mean = zoo::rollapply(delta_BAU_predict_mean,7,mean,align='right',fill=NA),
             d7_rollavg_delta_BAU_predict_sd = zoo::rollapply(delta_BAU_predict_sd,7,mean,align='right',fill=NA),
             d7_rollavg_observed_mean = zoo::rollapply(observed_mean,7,mean,align='right',fill=NA),
             d7_rollavg_BAU_mean = zoo::rollapply(BAU_mean,7,mean,align='right',fill=NA),
             d7_rollavg_CI_lower = zoo::rollapply(CI_lower,7,mean,align='right',fill=NA),
             d7_rollavg_CI_upper = zoo::rollapply(CI_upper,7,mean,align='right',fill=NA))
    
    return(predict_df_combo)
  }
}

#----------outputs model statistics in an easy to read format ---------------

urban_model_statistics = function(df, UK_code){
  names(df) = UK_code
  
  predict_df = plyr::ldply(df, 
                           data.frame) %>%
    filter(!is.na(date)) %>%
    rename(sites = .id)%>%
    select(sites, date, training_rsquared, training_mse, testing_rsquared, testing_mse) %>%
    mutate(training_rsquared = training_rsquared*100, testing_rsquared = testing_rsquared*100) %>%
    group_by(sites) %>%
    summarise_all(mean) %>%
    select(-date) %>%
    rename(Training_Rsquared=training_rsquared,
           Training_MSE=training_mse,
           Testing_Rsquared=testing_rsquared,
           Testing_MSE=testing_mse)%>%
    arrange(Testing_Rsquared) %>%
    ungroup() %>%
    mutate(site_id = row_number())%>%
    pivot_longer(-c(sites, site_id), names_to = "statistics") %>%
    separate(statistics, into = c("Data", "stat"), sep = "_")
  
  return(predict_df)
  
}
