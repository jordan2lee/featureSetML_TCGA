color_ht_col_labels <- function(membership_vector){
  #' Color of feature font text for heatmap columns names
  membership_vector[membership_vector == "" ] <- 'black'
  membership_vector[membership_vector == "+" ] <- 'red'
  return(membership_vector)
}
