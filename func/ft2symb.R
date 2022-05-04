ft2symb <- function(names_vector, platform){
  #' Convert GDAN-TMP ft name to gene symbol
  #' Each data platform has unique naming system
  library(dplyr)
  
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
      # clean up trailing ':' in names
      final_names <- c()
      for (g in name){
        if ( endsWith(g, ':)') ){
          final_names <- c(final_names, gsub(':','', g)  )
        } else {
          final_names <- c(final_names, g  )
        }
      }
    return(final_names)
  }
}
