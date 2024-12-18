determine_thresholds <- function(x, features, test_values, true_count){
  values <- function(y){
    def_features(x, features = features > y)[["metadata"]][["feature_id"]] |>
      unique() |>
      length() - 1
  }  
  
  plots <- ggplot()+
    geom_line(aes(x = test_values, y = values)) +
    geom_hline(aes(yintercept = true_count)) +
    #scale_x_log10() +
    scale_y_log10() +
    theme_bw()
  
  error <- values - true_count
  
  ideal_threshold <- test_values[which.min(abs(error))]
  
  ideal_threshold_error <- error[which.min(abs(error))]
  
  list(values, plots, ideal_threshold, ideal_threshold_error)
}
