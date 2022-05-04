full_cohort_name <- function(cancer){
  #' Cancer display names
  #' TCGA abbrev to full cancer type
  l1 <- list(
    'BRCA' = 'breast invasive carcinoma',
    'LGGGBM' = 'brain lower grade glioma and glioblastoma multiforme',
    'COADREAD' = 'colon adenocarcinoma and rectum adenocarcinoma',
    'SKCM' = 'skin cutaneous melanoma',
    'ACC' = 'adrenocortical carcinoma',
    'BLCA' = 'bladder urothelial carcinoma',
    'CESC' = 'cervical squamous cell carcinoma and endocervical adenocarcinoma',
    'ESCC' = 'esophageal squamous',
    'GEA' = 'gastroesophageal adenocarcinoma',
    'HNSC' = 'head and neck squamous cell carcinoma',
    'KIRCKICH' = 'kidney chromophobe and kidney renal clear cell carcinoma',
    'KIRP' = 'kidney renal papillary cell carcinoma',
    'LIHCCHOL' = 'liver hepatocellular carcinoma and cholangiocarcinoma',
    'LUAD' = 'lung adenocarcinoma',
    'LUSC' = 'lung squamous cell carcinoma',
    'MESO' = 'mesothelioma',
    'OV' = 'ovarian serous cystadenocarcinoma',
    'PAAD' = 'pancreatic adenocarcinoma',
    'PCPG' = 'pheochromocytoma and paraganglioma',
    'PRAD' = 'prostate adenocarcinoma',
    'SARC' = 'sarcoma',
    'TGCT' = 'testicular germ cell tumors',
    'THCA' = 'thyroid carcinoma',
    'THYM' = 'thymoma',
    'UCEC' = 'uterine corpus endometrial carcinoma',
    'UVM' = 'uveal melanoma'
  )
  if (is.null(l1[[cancer]]) == TRUE){
    return(cancer)
  }
  return(l1[[cancer]])
}
