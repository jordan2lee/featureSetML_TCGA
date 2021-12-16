full_cohort_name <- function(cancer){
  #' Display name of cancer for figures.
  l1 <- list(
    'BRCA' = 'Breast invasive carcinoma',
    'LGGGBM' = 'Brain lower grade glioma and glioblastoma multiforme',
    'COADREAD' = 'Colon adenocarcinoma and rectum adenocarcinoma',
    'SKCM' = 'Skin cutaneous melanoma',
    'ACC' = 'Adrenocortical carcinoma',
    'BLCA' = 'Bladder urothelial carcinoma',
    'CESC' = 'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
    'ESCC' = 'Esophageal squamous',
    'GEA' = 'Gastroesophageal adenocarcinoma',
    'HNSC' = 'Head and neck squamous cell carcinoma',
    'KIRCKICH' = 'Kidney chromophobe and kidney renal clear cell carcinoma',
    'KIRP' = 'Kidney renal papillary cell carcinoma',
    'LIHCCHOL' = 'Liver hepatocellular carcinoma and cholangiocarcinoma',
    'LUAD' = 'Lung adenocarcinoma',
    'LUSC' = 'Lung squamous cell carcinoma',
    'MESO' = 'Mesothelioma',
    'OV' = 'Ovarian serous cystadenocarcinoma',
    'PAAD' = 'Pancreatic adenocarcinoma',
    'PCPG' = 'Pheochromocytoma and paraganglioma',
    'PRAD' = 'Prostate adenocarcinoma',
    'SARC' = 'Sarcoma',
    'TGCT' = 'Testicular germ cell tumors',
    'THCA' = 'Thyroid carcinoma',
    'THYM' = 'Thymoma',
    'UCEC' = 'Uterine corpus endometrial carcinoma',
    'UVM' = 'Uveal melanoma'
  )
  if (is.null(l1[[cancer]]) == TRUE){
    return(cancer)
  }
  return(l1[[cancer]])
}
