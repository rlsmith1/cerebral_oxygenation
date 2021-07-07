
## A function to find the longest segment of the signal (where t >= 200s) with variance under 1000 ##

## Since variance does not strictly increase or decrease, the best way to find the longest list with a variance
## below a certain limit is to check all potential sublists, starting with the biggest ones and then shrinking
## until we find a smaller one. 

## our parameters: f_filter_signal(df, min_length = 200s*sampling rate, max_var = 1000)


f_filter_signal <- function(df, min_length, max_var){
  
  for (length in nrow(df):min_length) {
    
    for (start in 1:(nrow(df) - length + 1)) {
      
      df_temp <- df[start:(start + length - 1),]

      var <- df_temp$THC %>% var(na.rm = TRUE)
      
      if (var <= max_var) {
        
        return(df_temp)
        
      }
      
    }
    
  }
  
  return(NULL)
  
}




doParallel::registerDoParallel()
sampling_rate <- 1/.02

f1 <- f_filter_signal(df_raw_data %>% dplyr::filter(number == 1), min_length = 200*sampling_rate, max_var = 1000)




f_filter_signal <- function(df, min_length, max_var){
  
  for (length in nrow(df):min_length) {
    
    for (start in 1:(nrow(df) - length + 1)) {
      
      df_temp <- df[start:(start + length - 1),]

      var <- df_temp$roll_mean %>% var(na.rm = TRUE)
      
      if (var <= max_var) {
        
        return(df_temp)
        
      }
      
    }
    
  }
  
  return(NULL)
  
}



