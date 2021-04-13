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

get_colors_platform <- function(platform){
  #' Input the data platform. Example string "GEXP"
  #' Returns color code
  colors <- list(
    'MUTA' = '#00BFFF',
    'CNVR' = '#00688b',
    'METH' = '#43CD80',
    'GEXP' = '#FFA500',
    'MIR' = '#FF7F00'
  )
  return(colors[[platform]])
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
      platforms <- c('N:GEXP') # TODO del after dev
      # platforms <- c('N:METH', 'I:CNVR', 'B:MUTA', 'N:GEXP', 'N:MIR')
    } else if ( cancer == 'KIRCKICH' | cancer == 'LGGGBM' | cancer == 'LIHCCHOL'){
      platforms <- c('B:MUTA', 'N:METH', 'N:GEXP', 'I:CNVR')
    }
    return(platforms)
}

get_ymax_upset <- function(cancer){
  l1 <- list(
    'ACC' = 80, 'BLCA' = 70, 'BRCA' = 50,
    'CESC'= 60, 'COADREAD' = 70, 'ESCC' = 70,
    'GEA' = 70, 'HNSC' = 60, 'KIRP' = 70,
    'LGGGBM' = 70, 'LUAD' = 70, 'LUSC' = 60,
    'MESO' = 90, 'OV' = 80, 'PAAD' = 90,
    'PCPG' = 60, 'PRAD' = 60, 'SARC' = 50,
    'SKCM' = 70, 'TGCT' = 110, 'THCA' = 90,
    'THYM' = 90, 'UVM' = 100
  )
  return(l1[[cancer]])
}

title_info <- function(plat){
  l1 <- list(
    'GEXP' = 'Gene Expression',
    'MIR' = 'miRNA',
    'METH' = 'Methylation',
    'CNVR' = 'Copy Number Variation',
    'MUTA' = 'Mutation Status'
  )
  return(l1[[plat]])
}

get_gpar <- function(key){
  #' Input key and return value
  l1 <- list(
    'c' = 'black', # color
    'font_fam' = 'sans',
    'font_fam_ggplot' = 'Times New Roman', #'Arial',
    'main_title_size' = 18,
    'axis_size' = 16,
    'minor_axis_size' = 5,
    'legend_size_title'= 11,
    'legend_size'= 9,
    'annot_size' = 12, # ht bottom annot
    'model_overlap_size' = 9, # axis ticks
    'symbol_size' = 9,
    'pam_size' = 11
  )
  return(l1[[key]])
}
