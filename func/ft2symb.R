ft2symb <- function(names_vector, platform){
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
  } else if (platform == 'MUTA'){
      display_genes <- sapply(strsplit(names_vector, ':'), `[`,4)
      display_gene_details <- sapply(strsplit(names_vector, ':'), `[`,5)
      display_gene_type <- sapply(strsplit(names_vector, ':'), `[`,3)
      name <- paste(display_genes, ' (',display_gene_type, ':', display_gene_details, ')', sep='')
    return(name)
  }
}
