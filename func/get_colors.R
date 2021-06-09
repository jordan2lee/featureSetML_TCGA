get_colors <- function(df){
  #' Pull color codes based on total number of subtypes found in "Labels" column
  #' Input df that contains column "Labels" where subtypes are
  #' Returns vector of set heatmap color codes for subtypes
  nsubs <- df %>%
    select(Labels) %>%
    unique() %>%
    nrow() %>%
    as.integer()
  if (nsubs == 2){
    color_codes <- c(
      "1"='pink4',
      "2"='darkorange2'
    )
  } else if (nsubs == 3){
    color_codes <- c(
      "1"='pink4',
      "2"='darkorange2',
      "3"= 'sandybrown'
    )
  } else if (nsubs == 4){
    color_codes <- c(
      "1"='pink4',
      "2"='darkorange2',
      "3"= 'sandybrown',
      "4" = 'orange3'
    )
  } else if (nsubs == 5){
    color_codes <- c(
      "1"='tan',
      "2"='pink4',
      "3"='darkorange2',
      "4"= 'sandybrown',
      "5" = 'orange3'
    )
  } else if (nsubs == 6){
    color_codes <- c(
      "1" = 'saddlebrown',
      "2"='tan',
      "3"='pink4',
      "4"='darkorange2',
      "5"= 'sandybrown',
      "6" = 'orange3'
    )
  } else if (nsubs == 7){
    color_codes <- c(
      "1" = 'saddlebrown',
      "2" ='lightsalmon3',
      "3"='tan',
      "4"='pink4',
      "5"='darkorange2',
      "6"= 'sandybrown',
      "7" = 'orange3'
    )
  }
  return(color_codes)
}
