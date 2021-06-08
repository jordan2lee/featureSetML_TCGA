color_pam <- function(membership_vector){
  #' Color of feature column text to show PAM membership
  membership_vector[membership_vector == "" ] <- 'black'
  membership_vector[membership_vector == "+" ] <- 'red'
  return(membership_vector)
}

ft2gene_gexp <- function(names_vector, platform){
  library(dplyr)
  # Convert full ft name to gene symbol ONLY
  names_vector <- names_vector %>% as.vector()
  if (platform == 'GEXP'){
    display_genes <- sapply(strsplit(names_vector, '::'), `[`,2)
    display_genes <- sapply(strsplit(display_genes, ':'), `[`,1)
    return(display_genes)
  } else if (platform == 'METH'){
      display_genes <- sapply(strsplit(names_vector, ':'), `[`,4)
    return(display_genes)
  }
}
