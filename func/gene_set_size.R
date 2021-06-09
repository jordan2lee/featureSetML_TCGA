gene_set_size <- function(input_hallmark){
  # subset for hallmark rows
  tab<- mappings[mappings['gs_name']==input_hallmark,]
  # count unique gene symbols
  genes <- unique(unlist(tab['human_gene_symbol']))
  return(length(genes))
}
