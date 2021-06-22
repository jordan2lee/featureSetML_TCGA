jadbio_ft_file <- function(cancer){
  #' Input cancer and output the file containing jadbio ft importances
  #' TODO incorp this into the normal pipeline this is TEMP
  jadbio_files <- list(
    'LGGGBM' = 'src/jadbio_ft_importances_f1/LGGGBM_MULTIDATATYPE_markers_outcome_association_multisignature.csv',
    'BRCA' = 'src/jadbio_ft_importances_f1/BRCA_GEXP_markers_outcome_association.csv',
    'COADREAD' = 'src/jadbio_ft_importances_f1/COADREAD_MULTIDATATYPE_markers_outcome_association_multisignature.csv',
    'SKCM' = 'src/jadbio_ft_importances_f1/SKCM_MUTA_markers_outcome_association.csv'
  )
  return(jadbio_files[[cancer]])
}
