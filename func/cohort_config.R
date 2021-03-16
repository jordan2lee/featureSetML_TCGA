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


# get_platform_of_interest <- function(cancer){
#     #' Define platform to use for hallmark heatmap
#     #' Typically platform that subtypes are largely defined by
#     if (cancer == 'BRCA'){
#       platform_of_interest <- 'N:GEXP'
#     } else if (cancer == 'LGGGBM'){
#       platform_of_interest <- 'N:METH'
#     } else if (cancer == 'GEA'){
#       platform_of_interest <- 'N:METH'
#     }
#     return(platform_of_interest)
# }


get_platforms_present <- function(cancer){
    # All 5 platforms
    if (
      cancer == 'ACC' |
      cancer == 'BLCA' |
      cancer == 'BRCA' |
      cancer == 'CESC' |
      cancer == 'COADREAD' |
      cancer == 'ESCC' |
      cancer == 'GEA' |
      cancer == 'HNSC' |
      cancer == 'KIRP'|
      cancer == 'LUAD' |
      cancer == 'LUSC' |
      cancer == 'MESO' |
      cancer == 'OV' |
      cancer == 'PAAD' |
      cancer == 'PCPG' |
      cancer == 'PRAD' |
      cancer == 'SARC' |
      cancer == 'SKCM' |
      cancer == 'TGCT' |
      cancer == 'THCA' |
      cancer == 'THYM' |
      cancer == 'UCEC' |
      cancer == 'UVM' ){
      platforms <- c('N:METH', 'N:MIR', 'I:CNVR', 'B:MUTA', 'N:GEXP')
    } else if ( cancer == 'KIRCKICH' | cancer == 'LGGGBM' | cancer == 'LIHCCHOL'){
      platforms <- c('B:MUTA', 'N:METH', 'N:GEXP', 'I:CNVR')
    }
    return(platforms)
}
