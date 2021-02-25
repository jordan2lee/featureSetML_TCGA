#!/usr/bin/env Rscript

# install.packages('MXM')
# install.packages('rowr', quiet=TRUE)
# install.packages('argparse', quiet=TRUE)
suppressPackageStartupMessages(library(MXM))
suppressPackageStartupMessages(library(rowr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))

# create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--file", type='character', help="input file name with path")
parser$add_argument("-c", "--cancer", type='character', help='cancer such as brca, will be part of output naming')
parser$add_argument("-t", "--conditional_test", type='character', help='fbed test')
parser$add_argument("-m", "--method", type='character', help='fbed crit')
parser$add_argument("-k", "--kmax", type='integer', help='fbed k')
parser$add_argument("-th", "--thres", type='double', help='fbed thres')
parser$add_argument("-back", "--backphase", type='logical', default=TRUE, help='fbed backward phase, default is false')
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

######################################
# SUBSET FOR EACH PLATFORM (E.G. GEXP, MUTA, ETC)
#      CREATE LABEL AND FEATURE MATRICES  per platform
#####################################
names <- colnames(raw)

# Subset by PLATFORM = CREATE FEATURE MATRICES
ALL_ft_matrix <- raw

# CREATE LABEL MATRIX
#same for all platform subsets because didn't alter rows
label_matrix <- raw$Labels
label_matrix<- as.factor(label_matrix) # need as factor for FBED and SES

# Grab ft list
ALL_ft_matrix <- raw %>% select(-c(args$cancer, Labels))

######################################
# FEATURE SELECTION : ALL
#####################################
print('### starting ALL: fbed ebic ###')
setwd(args$output_path)


# Run algorithm
###########################
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

    ALL_ebic <- MXM::fbed.reg(target = label_matrix, dataset = ALL_ft_matrix, K = args$kmax, threshold = args$thres, test = args$conditional_test, method = args$method, backward = args$backphase)
    # Results + Output
    # 1. Create a list of selected features   -downstream summary .tsv
    ALL_select  <- c(colnames(ALL_ft_matrix[ALL_ebic$res[,'Vars']])) # fix this
    # print("selected ALL features")
    # print(ALL_select)

    #naMIng convention
    name = paste(args$cancer, "_", FtS, args$method, "--", platform, sep="")


    # 3. Summary Txt
    txtname = paste("summary-",name, "--combined_platform.txt", sep="")
    sink(txtname) #all info to be written below
    print(Sys.Date())
    print("##### ALL with FBED eBIC #####")
    print(paste('Cancer:', args$cancer, sep = ' '))
    print('Feature selection algorithm: FBED')
    print(paste('Method:', args$method, sep = ' '))
    print(paste('K:', args$kmax, sep = ' '))
    print(paste('Significance threshold:', args$thres, sep = ' '))
    print(paste('Test:', args$conditional_test, sep = ' '))
    print(paste('Backward selection phase:', args$backphase, sep = ' '))
    print("------------------------------------------------------")
    print("selected features")
    print(ALL_select)
    print("------------------------------------------------------")
    print("Selected <method> results:")
    print("R object: ALL_ebic$res")
    print('Cols = [Selected ft index, method, ft importance]')
    print(ALL_ebic$res)
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
    print("------------------------------------------------------")
    print('Runtime')
    print(ALL_ebic$runtime)
    sink()
}
print('saving file: ')
print(args$output_path)
print(txtname)

# 2. Create selected features tsv
# all <- cbind.fill(MUT_select, CNVR_select, GEXP_select, MI_select, MET_select, fill=NA)
# colnames(all) <- c("MUTA", "CNVR", "GEXP", "MI", "MET")
outputname = paste(args$cancer, "_", FtS, args$method, "--combined_platform.tsv", sep="")
write.table(ALL_select, file=outputname, quote=FALSE, sep=' ', col.names = TRUE, row.names=TRUE)
print('saving file: ')
print(args$output_path)
print(outputname)

print("SCRIPT COMPLETED FOR")
print(args$cancer)
