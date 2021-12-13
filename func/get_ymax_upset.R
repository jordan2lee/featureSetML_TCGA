get_ymax_upset <- function(cancer){
  #' Max y value for upset plot
  l1 <- list(
    'ACC' = 110, 'BLCA' = 70, 'BRCA' = 70,
    'CESC'= 60, 'COADREAD' = 50, 'ESCC' = 80,
    'GEA' = 50, 'HNSC' = 60, 'KIRCKICH' = 70, 'KIRP' = 80,
    'LGGGBM' = 60, 'LIHCCHOL' = 80, 'LUAD' = 80, 'LUSC' = 60,
    'MESO' = 90, 'OV' = 70, 'PAAD' = 30,
    'PCPG' = 60, 'PRAD' = 80, 'SARC' = 70,
    'SKCM' = 11.5, 'TGCT' = 70, 'THCA' = 70,
    'THYM' = 90, 'UVM' = 90
  )
  return(l1[[cancer]])
}
