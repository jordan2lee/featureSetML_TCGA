color_pam <- function(membership_vector){
  #' Color of feature column text to show PAM membership
  membership_vector[membership_vector == "" ] <- 'black'
  membership_vector[membership_vector == "+" ] <- 'red'
  return(membership_vector)
}
