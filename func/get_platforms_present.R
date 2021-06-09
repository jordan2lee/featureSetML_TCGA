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
