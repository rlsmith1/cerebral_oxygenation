
## A function to find the longest segment of the signal (where t >= 200s) within 0.5 standard deviations ##
## Purpose is to filter blips out of signal

## find the longest list with a standard deviation below a certain limit is to check all potential sublists

## our parameters: f_filter_signal(df, max_sd = 0.5)


f_filter_signal <- function(df, max_sd) {
  
  start = 1
  end = 1
  max_start = 1
  max_end = 0
  
  df <- df %>% filter(!is.na(roll_mean))
  
  mean <- df$roll_mean %>% mean(na.rm = TRUE)
  sd <- df$roll_mean %>% sd(na.rm = TRUE)
  
  for (i in 1:nrow(df)) {
    
    if ((mean - max_sd*sd) < df[i,]$roll_mean && df[i,]$roll_mean < (mean + max_sd*sd)) {
      
      end <- i
      
    } else {
      
      start <- i
      
    }
    
    if ((end - start) > (max_end - max_start)) {
      
      max_end <- end
      max_start <- start
      
    }
    
  }
  
  return(df[max_start:max_end,])
  
}





## Alternative function where you can filter with a max variance

f_filter_signal_var <- function(df, min_length, max_var){
  
  for (length in nrow(df):min_length) {
    
    for (start in 1:(nrow(df) - length + 1)) {
      
      df_temp <- df[start:(start + length - 1),]
      df_temp <- df_temp %>% mutate(start = start, length = length)
      
      var <- df_temp$roll_mean %>% var(na.rm = TRUE)
      
      if (var <= max_var) {
        
        return(df_temp)
        
      }
      
    }
    
  }
  
  return(NULL)
  
}








