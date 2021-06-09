build_hallmark_vect <- function(hallmark, ftnames_order, platform_of_interest){
  #' Create hallmark vector
  #' for builidng hallmark heatmap
  fts_checked <- c() # sanity check
  hallmark_present <- c()
  for (feature in ftnames_order){
    # 1. Preprocess - to gene symbol
    if (platform_of_interest == 'N:GEXP'){
      GENE <- unlist(strsplit(feature, '::'))[2]
      GENE <- unlist(strsplit(GENE, ':'))[1]
  } else if (platform_of_interest == 'N:METH' || platform_of_interest == 'B:MUTA' || platform_of_interest == 'I:CNVR'){
      GENE <- unlist(strsplit(feature, ':'))[4]
    }
    # 2. Hallmark Mapping
    halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>%
      as.vector()
    # If hallmark present
    if (hallmark %in% halls){
      i <- match(hallmark, halls)
      hallmark_present <- c(hallmark_present, 1)
      fts_checked <- c(fts_checked, feature)
    }
    else {
      hallmark_present <- c(hallmark_present, 0)
      fts_checked <- c(fts_checked, feature)
    }
  }
  return(hallmark_present)
}
