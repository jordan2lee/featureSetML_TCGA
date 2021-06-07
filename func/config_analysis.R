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

get_platforms_present <- function(cancer){
  #' Specifies which platforms present for each cancer type
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
    platforms <- c('N:METH', 'I:CNVR', 'B:MUTA', 'N:GEXP', 'N:MIR')
  } else if ( cancer == 'KIRCKICH' | cancer == 'LGGGBM' | cancer == 'LIHCCHOL'){
    platforms <- c('B:MUTA', 'N:METH', 'N:GEXP', 'I:CNVR')
  }
  return(platforms)
}

get_ymax_upset <- function(cancer){
  #' Max y value for upset plot
  l1 <- list(
    'ACC' = 110, 'BLCA' = 70, 'BRCA' = 70,
    'CESC'= 60, 'COADREAD' = 50, 'ESCC' = 80,
    'GEA' = 50, 'HNSC' = 60, 'KIRP' = 80,
    'LGGGBM' = 60, 'LUAD' = 80, 'LUSC' = 60,
    'MESO' = 90, 'OV' = 70, 'PAAD' = 30,
    'PCPG' = 60, 'PRAD' = 80, 'SARC' = 70,
    'SKCM' = 11.5, 'TGCT' = 70, 'THCA' = 70,
    'THYM' = 90, 'UVM' = 90
  )
  return(l1[[cancer]])
}

title_info <- function(plat){
  #' Display name of platform for figures
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
    'c' = 'black', # font color
    'font_fam' = 'sans',
    'font_fam_ggplot' = 'Arial',
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

platform_display_text <- function(platform){
  #input platform and convert to the text for figure display
  l1 <- list(
    'GEXP' = 'mRNA',
    'CNVR' = 'CN',
    'METH' = 'DNA Methylation',
    'MIR' = 'miR'
  )
  return(l1[[platform]])
}

ht_matrix_colors <- function(platform){
  l1 <- list(
    'CNVR' = structure(c('blue', 'white', 'red'), names = c(-1, 0, 1)),
    'MIR' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'METH' = NULL,
    'GEXP' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'MUTA' = structure(c('blue','red'), names = c(0, 1))
  )
  return(l1[[platform]])
}
