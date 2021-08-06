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
      "1" = '#9D8562',
      "2" ='#97ABC6'
    )
  } else if (nsubs == 3){
    color_codes <- c(
      "1" = '#9D8562',
      "2" ='#97ABC6',
      "3"='#97C6B2'
    )
  } else if (nsubs == 4){
    color_codes <- c(
      "1" = '#9D8562',
      "2" ='#97ABC6',
      "3"='#97C6B2',
      "4"='#E4E08D'
    )
  } else if (nsubs == 5){
    color_codes <- c(
      "1" = '#9D8562',
      "2" ='#97ABC6',
      "3"='#97C6B2',
      "4"='#E4E08D',
      "5"='#C5A8AB'
    )
  } else if (nsubs == 6){
    color_codes <- c(
      "1" = '#9D8562',
      "2" ='#97ABC6',
      "3"='#97C6B2',
      "4"='#E4E08D',
      "5"='#C5A8AB',
      "6"= '#C9D8E4'
    )
  } else if (nsubs == 7){
    color_codes <- c(
      "1" = '#9D8562',
      "2" ='#97ABC6',
      "3"='#97C6B2',
      "4"='#E4E08D',
      "5"='#C5A8AB',
      "6"= '#C9D8E4',
      "7" = '#878791'
    )
  }
  return(color_codes)
}
