subtype_common_name <- function(cancer, subtype_matrix, conversion_vector){
  #' Display name of cancer subtype for figures.
  l1 <- list(
    'BRCA' = c(
      '1' = 'LumA',
      '2' = 'LumB',
      '3' = 'Basal',
      '4' = 'Her2'
    ),
    'COADREAD' = c(
      '1' = 'CIN and non-CIMP',
      '2' = 'CIN and CIMP low',
      '3' = 'GS and non-CIMP high',
      '4' = 'MSI and CIMP High',
      '5' = 'MSI and non-CIMP High',
      '6' = 'MSS and CIMP high'#,
      # '7' = 'HM SNV (POLE)' # TODO unclear if this subtype was excluded from my input data
    ),
    'LGGGBM' = c(
      '1' = 'G CIMP high',
      '2' = 'Mesenchymal like',
      '3' = 'Codel',
      '4' = 'Classic like',
      '5' = 'LGm6 GBM',
      '6' = 'PA like',
      '7' = 'G CIMP low'
    ),
    'SKCM' = c(
      '1' = 'BRAF hotspot',
      '2' = 'RAS hotspot',
      '3' = 'Triple WT',
      '4' = 'NF1 LoF'
    )
  )
  if (conversion_vector == 'yes'){
    return(l1[cancer] %>% as.data.frame())
  } else {
    m = as.vector(
      lapply(
        subtype_matrix,
        FUN = function(x) l1[[cancer]][x] %>% as.list() %>% unlist() %>% unname()
      )
    )
    return(unlist(m))
  }
}
