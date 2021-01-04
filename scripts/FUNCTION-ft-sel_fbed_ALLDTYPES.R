#!/usr/bin/env Rscript

# install.packages('MXM')
# install.packages('rowr', quiet=TRUE)
# install.packages('argparse', quiet=TRUE)
suppressPackageStartupMessages(library(MXM))
suppressPackageStartupMessages(library(rowr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))


#### hardcoded
inputdate = 'ran on 12/9/20'
####

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

# # Create list of col headers
# MUT_names <- names[grepl(":MUTA:", names)]  #list of feature headers
# CNVR_names <- names[grepl(":CNVR::", names)]  #list of feature headers
# GEXP_names <- names[grepl(":GEXP::", names)]  #list of feature headers
# MI_names <- names[grepl(":MIR::", names)]  #list of feature headers
# MET_names <- names[grepl(":METH:", names)]  #list of feature headers
# print("sanity check - these should be equal: ")
# print(length(MUT_names)+length(CNVR_names)+length(GEXP_names)+length(MI_names)+length(MET_names)+2)
# print(dim(raw)[2])

# Subset by PLATFORM = CREATE FEATURE MATRICES
ALL_ft_matrix <- raw
# MUT_ft_matrix <- raw[MUT_names] #table were all fts from only one platform
# CNVR_ft_matrix <- raw[CNVR_names] #table were all fts from only one platform
# GEXP_ft_matrix <- raw[GEXP_names] #table were all fts from only one platform
# MI_ft_matrix <- raw[MI_names] #table were all fts from only one platform
# MET_ft_matrix <- raw[MET_names] #table were all fts from only one platform

# CREATE LABEL MATRIX
#same for all platform subsets because didn't alter rows
label_matrix <- raw$Labels
label_matrix<- as.factor(label_matrix) # need as factor for FBED and SES

# Grab ft list
ALL_ft_matrix <- raw %>% select(-c(args$class, Labels))
######################################
# FEATURE SELECTION : ALL
#####################################
print('### starting ALL: fbed ebic ###')
setwd(args$output_path)


# Run algorithm
###########################
class = args$class
thres = 0.05
conditional_test = "testIndMultinom"
kmax = 2
crit = 'eBIC'
FtS = "fbed"
platform = 'ALL'
###########################
#If not features of this platform of raw data, then DON'T run feature selection
if (dim(ALL_ft_matrix)[2] == 0){
    print("no input features of this platform, therefore skipping")
    ALL_select <- c()
}
#If features then run feature selection
if (dim(ALL_ft_matrix)[2] != 0){

    ALL_ebic <- MXM::fbed.reg(target = label_matrix, dataset = ALL_ft_matrix, K = kmax, threshold = thres, test = conditional_test, method = crit, backward = TRUE)
    # print("algorithm done and now writing results")
    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    ALL_select  <- c(colnames(ALL_ft_matrix[ALL_ebic$res[,'Vars']])) # fix this
    # print("selected ALL features")
    # print(ALL_select)

    #naMIng convention
    name = paste(class, "_", FtS, crit, "--", platform, sep="")


    # 3. Summary Txt
    txtname = paste("summary-",name, "--combined_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(inputdate)
    print(args$class)
    print("##### ALL with FBED eBIC #####")
    print("hyperparameters used (max k, threshold, test name)")
    print(class)
    print(crit)
    print(kmax)
    print(thres)
    print(conditional_test)
    print("------------------------------------------------------")
    print(ALL_ebic$runtime)
    print("------------------------------------------------------")
    print("eBIC results:")
    print("ALL_ebic$res")
    print(ALL_ebic$res)
    print("------------------------------------------------------")
    print("selected features")
    print(ALL_select)
    print("------------------------------------------------------")
    print("selected feature scores")
    print(ALL_ebic$res[,2])
    print("------------------------------------------------------")
    print("info matrix: FORWARD PHASE")
    print(ALL_ebic$info)
    print("------------------------------------------------------")
    print("info matrix: BACKWARD PHASE")
    print("shows variables removed from bkward phase")
    print(ALL_ebic$back.rem)
    print("------------------------------------------------------")
    print("number of models that were fitted to the backward phase")
    print(ALL_ebic$back.n.tests)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

# 2. Create selected features tsv
# all <- cbind.fill(MUT_select, CNVR_select, GEXP_select, MI_select, MET_select, fill=NA)
# colnames(all) <- c("MUTA", "CNVR", "GEXP", "MI", "MET")
outputname = paste(class, "_", FtS, crit, "--combined_platform.tsv", sep="")
write.table(ALL_select, file=outputname, quote=FALSE, sep='\t', col.names = TRUE, row.names=TRUE)
print('saving file: ')
print(args$output_path)
print(outputname)

print("SCRIPT COMPLETED FOR")
print(args$class)
