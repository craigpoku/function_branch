



#-----------creates linear function based on test data setup by the user--------------------

linear_function_generator = function(index, df, seq, epsilon){
  f = approxfun(x = index, y = df)
  seq_range_x = c(1:seq) 
  linear_fit = f(seq_range_x)
  
  f_x_df = data.frame(seq_range_x, linear_fit) %>%
    rename(index = seq_range_x, value = linear_fit) %>%
    mutate(value = jitter(value, factor=epsilon, amount = NULL))%>%
    drop_na()
  
  return(f_x_df)
}

#-----------Updated function that determines the CP based on ratio of previous to current point using the 2nd derivative of a univariant time series ----------

multi_var_ts_gradient_cp_detection = function(df, window_length_vector, df_header_code, 
                                              epsilon, cp_factor, date = TRUE){
  
  df_filter = df %>%
    filter(df_header == df_header_code)
  
  if(date==TRUE){
    roll_regression = rollRegres::roll_regres(value ~ date, df_filter, 
                                              width = window_length_vector,
                                              do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))  
    
    
    roll_reformat_cp = as.data.frame(roll_regression$coefs) %>%
      rename(grad = date) %>%
      mutate(date = df_filter$date,
             r.squareds = roll_regression$r.squareds,
             data = df_filter$value,
             df_label = df_header_code,
             window_length_level = as.factor(window_length_vector),
             derv_2nd = as.numeric(abs(pracma::gradient(grad))),
             cp = derv_2nd/lag(derv_2nd) > cp_factor &
               lag(derv_2nd) > epsilon,
      ) %>%rename("Input dataset" = data,
                  "Rolling gradient" = grad,
                  "2nd derivative" = derv_2nd)%>%
      select(-"(Intercept)") %>%
      drop_na() %>%
      pivot_longer(-c(date, window_length_level,df_label, cp), 
                   names_to = "variables")%>%
      mutate(variables = factor(variables, 
                                levels = c("Input dataset", "Rolling gradient", "2nd derivative"
                                           , "r.squareds")))
  }
  else{
    roll_regression = rollRegres::roll_regres(value ~ index, df_filter, 
                                              width = window_length_vector,
                                              do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))  
    
    
    roll_reformat_cp = as.data.frame(roll_regression$coefs) %>%
      rename(grad = index) %>%
      mutate(index = df_filter$index,
             r.squareds = roll_regression$r.squareds,
             data = df_filter$value,
             df_label = df_header_code,
             window_length_level = as.factor(window_length_vector),
             derv_2nd = as.numeric(abs(pracma::gradient(grad))),
             cp = derv_2nd/lag(derv_2nd) > cp_factor &
               derv_2nd-lag(derv_2nd) > 0 &
               lag(derv_2nd) > epsilon
      ) %>%rename("Test dataset" = data,
                  "Rolling gradient" = grad,
                  "2nd derivative" = derv_2nd)%>%
      select(-"(Intercept)") %>%
      drop_na() %>%
      pivot_longer(-c(index, window_length_level,df_label, cp), 
                   names_to = "variables")%>%
      mutate(variables = factor(variables, 
                                levels = c("Test dataset", "Rolling gradient", "2nd derivative"
                                           , "r.squareds")))
    
    
  }
  return(roll_reformat_cp)
}


#--------- function to test what is a suitable choice in window length --------

change_point_model_statistics = function(df, window_length, stats = TRUE){
  df_reformat = df %>%
    filter(variables == "Test dataset", window_length_level == window_length) %>%
    mutate(date = lubridate::as_date(date))
  
  cp_df  = df_reformat %>%
    filter(cp == TRUE)
  
  ApproxFun = approxfun(x = cp_df$date, y = cp_df$value)
  Dates_seq <- seq.Date(ymd(df$date[1]), ymd(tail(df$date, n=1)), by = 1)
  LinearFit = ApproxFun(Dates_seq)
  
  df_approx = data.frame(Dates_seq, LinearFit)%>%
    rename(date = Dates_seq, approx_value = LinearFit)
  
  df_new = df %>%
    filter(window_length_level == window_length, variables == "Test dataset") 
  
  if(stats == TRUE){
    
    df_combo_temp = left_join(df_new, df_approx, by = "date") %>%
      drop_na() %>%
      select(date, value, approx_value)
    
    rmse = Metrics::rmse(df_combo_temp$value, df_combo_temp$approx_value)
    rsquared = cor(df_combo_temp$value, df_combo_temp$approx_value, method ="pearson")
    df_combo = data.frame(window_length, rmse, rsquared)
  }
  else{
    df_combo = left_join(df_new, df_approx, by = "date")%>%
      mutate(cp = df_reformat$cp)  %>%
      select(date, window_length_level, cp, value, approx_value)  %>%
      drop_na() %>%
      rename("f(x)" = value, "Approxed f(x)" = approx_value) %>%
      pivot_longer(-c(date, window_length_level, cp), names_to = "variables")
  }
  
  return(df_combo)
}