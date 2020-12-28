#!/usr/bin/env Rscript
# install.packages('MXM')
library(MXM)
# install.packages('rowr', quiet=TRUE)
library(rowr)
# install.packages('argparse', quiet=TRUE)
library(argparse)


# create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--file", type='character', help="input file name with path")
parser$add_argument("-c", "--class", type='character', help='class such as brca, will be part of output naming')
parser$add_argument('-op', '--output_path', type='character', help='output path')
args <- parser$parse_args()


######################################
# SETUP
#####################################
# set seed
set.seed(1)

# Read in Data
raw  = read.table(args$file, sep='\t', header=TRUE, check.names = FALSE)

# EDA: Inspect dataset
dim(raw)
# #Ensure no NA or NaN in dataset
# message("NA found is = ", sum(is.na(raw))) #check NA = 0
# message("NaN found is = ", sum(is.nan(as.matrix(raw)))) #chek NaN = 0

######################################
# SUBSET FOR EACH PLATFORM (E.G. GEXP, MUTA, ETC)
#      CREATE LABEL AND FEATURE MATRICES  per platform
#####################################
names <- colnames(raw)

# Create list of col headers
MUT_names <- names[grepl(":MUTA:", names)]  #list of feature headers
CNVR_names <- names[grepl(":CNVR::", names)]  #list of feature headers
GEXP_names <- names[grepl(":GEXP::", names)]  #list of feature headers
MI_names <- names[grepl(":MIR::", names)]  #list of feature headers
MET_names <- names[grepl(":METH:", names)]  #list of feature headers
print("sanity check - these should be equal: ")
print(length(MUT_names)+length(CNVR_names)+length(GEXP_names)+length(MI_names)+length(MET_names)+2)
print(dim(raw)[2])

# Subset by PLATFORM = CREATE FEATURE MATRICES
MUT_ft_matrix <- raw[MUT_names] #table were all fts from only one platform
CNVR_ft_matrix <- raw[CNVR_names] #table were all fts from only one platform
GEXP_ft_matrix <- raw[GEXP_names] #table were all fts from only one platform
MI_ft_matrix <- raw[MI_names] #table were all fts from only one platform
MET_ft_matrix <- raw[MET_names] #table were all fts from only one platform

# CREATE LABEL MATRIX
#same for all platform subsets because didn't alter rows
label_matrix <- raw$Labels
label_matrix<- as.factor(label_matrix) # need as factor for FBED and SES

######################################
# FEATURE SELECTION : MUTA
#####################################
print('### starting MUTA: fbed ebic ###')
setwd(args$output_path)


# Run algorithm
###########################
class = args$class
thres = 0.05
conditional_test = "testIndMultinom"
kmax = 2
crit = 'eBIC'
FtS = "fbed"
platform = 'MUT'
###########################
#If not features of this platform of raw data, then DON'T run feature selection
if (dim(MUT_ft_matrix)[2] == 0){
    print("no input features of this platform, therefore skipping")
    MUT_select <- c()
}
#If features then run feature selection
if (dim(MUT_ft_matrix)[2] != 0){

    MUT_ebic <- MXM::fbed.reg(target = label_matrix, dataset = MUT_ft_matrix, K = kmax, threshold = thres, test = conditional_test, method = crit, backward = TRUE)
    # print("algorithm done and now writing results")
    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    MUT_select  <- c(colnames(MUT_ft_matrix[MUT_ebic$res[,1]]))
    # print("selected MUT features")
    # print(MUT_select)

    #naMIng convention
    name = paste(class, "_", FtS, crit, "--", platform, sep="")


    # 3. Summary Txt
    txtname = paste("summary-",name, "--indep_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(args$class)
    print("##### MUTA with FBED eBIC #####")
    print("hyperparameters used (max k, threshold, test name)")
    print(class)
    print(crit)
    print(kmax)
    print(thres)
    print(conditional_test)
    print("------------------------------------------------------")
    print(MUT_ebic$runtime)
    print("------------------------------------------------------")
    print("eBIC results:")
    print("MUT_ebic$res")
    print(MUT_ebic$res)
    print("------------------------------------------------------")
    print("selected features")
    print(MUT_select)
    print("------------------------------------------------------")
    print("selected feature scores")
    print(MUT_ebic$res[,2])
    print("------------------------------------------------------")
    print("info matrix: FORWARD PHASE")
    print(MUT_ebic$info)
    print("------------------------------------------------------")
    print("info matrix: BACKWARD PHASE")
    print("shows variables removed from bkward phase")
    print(MUT_ebic$back.rem)
    print("------------------------------------------------------")
    print("number of models that were fitted to the backward phase")
    print(MUT_ebic$back.n.tests)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

######################################
# FEATURE SELECTION : CNVR
#####################################
print('### starting CNVR: fbed ebic ###')

# Run algorithm
###########################
class = args$class
thres = 0.05
conditional_test = "testIndMultinom"
kmax = 2
crit = 'eBIC'
platform = "CNVR"
###########################
#If not features of this platform of raw data, then DON'T run feature selection
if (dim(CNVR_ft_matrix)[2] == 0){
    print("no input features of this platform, therefore skipping")
    CNVR_select <- c()
}
#If features then run feature selection
if (dim(CNVR_ft_matrix)[2] != 0){

    CNVR_ebic <- MXM::fbed.reg(target = label_matrix, dataset = CNVR_ft_matrix, K = kmax, threshold = thres, test = conditional_test, method = crit, backward = TRUE)
    # print("algorithm done and now writing results")

    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    CNVR_select  <- c(colnames(CNVR_ft_matrix[CNVR_ebic$res[,1]]))
    # print("selected CNVR features")
    # print(CNVR_select)

    #naming convention
    name = paste(class, "_", FtS, crit, "--", platform, sep="")

    # 3. Summary Txt
    txtname = paste("summary-",name, "--indep_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(args$class)
    print("##### CNVR with FBED eBIC #####")
    print("hyperparameters used (max k, threshold, test name)")
    print(class)
    print(crit)
    print(kmax)
    print(thres)
    print(conditional_test)
    print("------------------------------------------------------")
    print(CNVR_ebic$runtime)
    print("------------------------------------------------------")
    print("eBIC results:")
    print("CNVR_ebic$res")
    print(CNVR_ebic$res)
    print("------------------------------------------------------")
    print("selected features")
    print(CNVR_select)
    print("------------------------------------------------------")
    print("selected feature scores")
    print(CNVR_ebic$res[,2])
    print("------------------------------------------------------")
    print("info matrix: FORWARD PHASE")
    print(CNVR_ebic$info)
    print("------------------------------------------------------")
    print("info matrix: BACKWARD PHASE")
    print("shows variables removed from bkward phase")
    print(CNVR_ebic$back.rem)
    print("------------------------------------------------------")
    print("number of models that were fitted to the backward phase")
    print(CNVR_ebic$back.n.tests)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

######################################
# FEATURE SELECTION : GEXP
#####################################
print('### starting GEXP: fbed ebic ###')

# Run algorithm
###########################
class = args$class
thres = 0.05
conditional_test = "testIndMultinom"
kmax = 2
crit = 'eBIC'
platform = "GEXP"
###########################

#If not features of this platform of raw data, then DON'T run feature selection
if (dim(GEXP_ft_matrix)[2] == 0){
    print("no input features of this platform, therefore skipping")
    GEXP_select <- c()
}
#If features then run feature selection
if (dim(GEXP_ft_matrix)[2] != 0){

    GEXP_ebic <- MXM::fbed.reg(target = label_matrix, dataset = GEXP_ft_matrix, K = kmax, threshold = thres, test = conditional_test, method = crit, backward = TRUE)
    # print("algorithm done and now writing results")

    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    GEXP_select  <- c(colnames(GEXP_ft_matrix[GEXP_ebic$res[,1]]))
    # print("selected GEXP features")
    # print(GEXP_select)

    #naming convention
    name = paste(class, "_", FtS, crit, "--", platform, sep="")

    # 3. Summary Txt
    txtname = paste("summary-",name, "--indep_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(args$class)
    print("##### GEXP with FBED eBIC #####")
    print("hyperparameters used (max k, threshold, test name)")
    print(class)
    print(crit)
    print(kmax)
    print(thres)
    print(conditional_test)
    print("------------------------------------------------------")
    print(GEXP_ebic$runtime)
    print("------------------------------------------------------")
    print("eBIC results:")
    print("GEXP_ebic$res")
    print(GEXP_ebic$res)
    print("------------------------------------------------------")
    print("selected features")
    print(GEXP_select)
    print("------------------------------------------------------")
    print("selected feature scores")
    print(GEXP_ebic$res[,2])
    print("------------------------------------------------------")
    print("info matrix: FORWARD PHASE")
    print(GEXP_ebic$info)
    print("------------------------------------------------------")
    print("info matrix: BACKWARD PHASE")
    print("shows variables removed from bkward phase")
    print(GEXP_ebic$back.rem)
    print("------------------------------------------------------")
    print("number of models that were fitted to the backward phase")
    print(GEXP_ebic$back.n.tests)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

######################################
# FEATURE SELECTION : MI
#####################################
print('### starting MI: fbed ebic ###')

# Run algorithm
###########################
class = args$class
thres = 0.05
conditional_test = "testIndMultinom"
kmax = 2
crit = 'eBIC'
platform = "MI"
###########################
#If not features of this platform of raw data, then DON'T run feature selection
if (dim(MI_ft_matrix)[2] == 0){
    print("no input features of this platform, therefore skipping")
    MI_select <- c()
}
#If features then run feature selection
if (dim(MI_ft_matrix)[2] != 0){

    MI_ebic <- MXM::fbed.reg(target = label_matrix, dataset = MI_ft_matrix, K = kmax, threshold = thres, test = conditional_test, method = crit, backward = TRUE)
    # print("algorithm done and now writing results")

    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    MI_select  <- c(colnames(MI_ft_matrix[MI_ebic$res[,1]]))
    # print("selected MI features")
    # print(MI_select)

    #naming convention
    name = paste(class, "_", FtS, crit, "--", platform, sep="")

    # 3. Summary Txt
    txtname = paste("summary-",name, "--indep_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(args$class)
    print("##### MI with FBED eBIC #####")
    print("hyperparameters used (max k, threshold, test name)")
    print(class)
    print(crit)
    print(kmax)
    print(thres)
    print(conditional_test)
    print("------------------------------------------------------")
    print(MI_ebic$runtime)
    print("------------------------------------------------------")
    print("eBIC results:")
    print("MI_ebic$res")
    print(MI_ebic$res)
    print("------------------------------------------------------")
    print("selected features")
    print(MI_select)
    print("------------------------------------------------------")
    print("selected feature scores")
    print(MI_ebic$res[,2])
    print("------------------------------------------------------")
    print("info matrix: FORWARD PHASE")
    print(MI_ebic$info)
    print("------------------------------------------------------")
    print("info matrix: BACKWARD PHASE")
    print("shows variables removed from bkward phase")
    print(MI_ebic$back.rem)
    print("------------------------------------------------------")
    print("number of models that were fitted to the backward phase")
    print(MI_ebic$back.n.tests)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

######################################
# FEATURE SELECTION : MET
#####################################
print('### starting MET: fbed ebic ###')

# Run algorithm
###########################
class = args$class
thres = 0.05
conditional_test = "testIndMultinom"
kmax = 2
crit = 'eBIC'
platform = "MET"
###########################
#If not features of this platform of raw data, then DON'T run feature selection
if (dim(MET_ft_matrix)[2] == 0){
    print("no input features of this platform, therefore skipping")
    MET_select <- c()
}
#If features then run feature selection
if (dim(MET_ft_matrix)[2] != 0){
    MET_ebic <- MXM::fbed.reg(target = label_matrix, dataset = MET_ft_matrix, K = kmax, threshold = thres, test = conditional_test, method = crit, backward = TRUE)
    # print("algorithm done and now writing results")

    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    MET_select  <- c(colnames(MET_ft_matrix[MET_ebic$res[,1]]))
    # print("selected MET features")
    # print(MET_select)

    #naming convention
    name = paste(class, "_", FtS, crit, "--", platform, sep="")

    # 3. Summary Txt
    txtname = paste("summary-",name, "--indep_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(args$class)
    print("##### MET with FBED eBIC #####")
    print("hyperparameters used (max k, threshold, test name)")
    print(class)
    print(crit)
    print(kmax)
    print(thres)
    print(conditional_test)
    print("------------------------------------------------------")
    print(MET_ebic$runtime)
    print("------------------------------------------------------")
    print("eBIC results:")
    print("MET_ebic$res")
    print(MET_ebic$res)
    print("------------------------------------------------------")
    print("selected features")
    print(MET_select)
    print("------------------------------------------------------")
    print("selected feature scores")
    print(MET_ebic$res[,2])
    print("------------------------------------------------------")
    print("info matrix: FORWARD PHASE")
    print(MET_ebic$info)
    print("------------------------------------------------------")
    print("info matrix: BACKWARD PHASE")
    print("shows variables removed from bkward phase")
    print(MET_ebic$back.rem)
    print("------------------------------------------------------")
    print("number of models that were fitted to the backward phase")
    print(MET_ebic$back.n.tests)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

# 2. Create selected features tsv
all <- cbind.fill(MUT_select, CNVR_select, GEXP_select, MI_select, MET_select, fill=NA)
colnames(all) <- c("MUTA", "CNVR", "GEXP", "MI", "MET")
#write.table(all, file='test.tsv', quote=FALSE, sep='\t', col.names = TRUE, row.names=TRUE)
outputname = paste(class, "_", FtS, crit, "--indep_platform.tsv", sep="")
write.table(all, file=outputname, quote=FALSE, sep='\t', col.names = TRUE, row.names=TRUE)
print('saving file: ')
print(args$output_path)
print(outputname)

print("SCRIPT COMPLETED FOR")
print(args$class)
