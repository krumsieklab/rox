# This scripts downloads, preprocesses, and saves the 6 metabolomics datasets 
# used in the `recovering high-confidence hits` analysis of the paper 
#
# The execution takes about 10-15min on an Apple M1 Pro

library(dplyr)
library(magrittr)
library(readxl)
# maplet
if (!("maplet" %in% rownames(installed.packages()))) stop("This script needs the maplet package. devtools::install_github(repo='krumsieklab/maplet', subdir='maplet')")
library(maplet)
 

# Initialize dataset lists
datasets <- list()

# CONVENTION:
# all datasets are identified by their list name (key) and must contain the following fields;
# $X        - data matrix
# $X.minimp - data matrix with missing values imputed using minimum value
# $X.knnimp - data matrix with missing values imputed using knn imputation
# $Covars   - Covariate matrix
# $groups   - two-group vector of sample classes (e.g. diabetes yes/no)


#### Functions ----

# generate maplet object
mapletify_data <- function(assay, colData, rowData){
  
  # generate Summarized Experiment
  SummarizedExperiment(assays = assay %>% t,
                       colData = colData,
                       rowData = rowData) %>% {.}
    # validate SE to use with maplet
    # maplet::mt_clean_validate_se()
}

# preprocessing steps
fun_preprocessing <- function(D){
  
  # normalize
  D %<>%
    # probabilistic quotient normalization
    maplet::mt_pre_norm_quot(feat_max = 0.2) %>%
    # log2-transform
    maplet:: mt_pre_trans_log()
  
  # knn imputation
  D_knn <- D %>%
    maplet::mt_pre_impute_knn()
  
  # minimum imputation
  D_min <- D %>%
    maplet::mt_pre_impute_min()
  
  # return all three SEs
  list(norm = D,
       knnimp = D_knn,
       minimp = D_min)
}

# load precalculated data from file
cached_block <- function(
  file,          # file name to save to or read from
  checksum=NA,   # file checksum to check
  verbose=F,     # print messages to console?
  block          # block of code {} 
) {
  
  # check if cached file exists
  if (!file.exists(file)) {
    # file does not exist
    # need to execute code and store result in file
    if (verbose) cat("Executing code block.\n")
    res <- block
    if (verbose) cat(sprintf("Saving results to '%s'.\n", file))
    save(res, file=file)
    # return results
    res
  } else {
    # file exists
    # verify checksum?
    if (!is.na(checksum)) {
      # calculate real checksum
      md5 <- tools::md5sum(file)[[1]]
      # crash if wrong
      if (checksum != md5) stop(sprintf("Wrong checksum for %s, expected: %s, actual: %s\n", file, checksum, md5))
    }
    # load precalculated data
    if (verbose) cat(sprintf("Loading results from '%s'.\n", file))
    load(file)
    # return
    res
    
  }
}

# Downloads files from the web and verifies their checksum. Will use local copy in current directory, if it exists
load.web.file <- function(
  url, md5sum, outfile, zipfile = F) {
  # check if local file exists
  if (file.exists(outfile)) {
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) stop(sprintf("Local file %s has wrong checksum: %s", outfile, realsum))
    # do not delete wrong file, it was already here before
    
  } else {
    if(zipfile){
      # download file
      temp <- tempfile()
      download.file(url,temp)
      unzip(zipfile = temp, files = outfile, exdir = ".")
    } else {
      # download file
      download.file(url, outfile)
    }
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) { 
      # delete wrong file
      unlink(outfile)
      stop(sprintf("Remote file %s has wrong checksum: %s", url, realsum))
    }
  }
}



#### Load QMdiab Data ----

# Original Publication: https://academic.oup.com/bioinformatics/article/35/3/532/5056040#supplementary-data

# cached block, to avoid recomputing of KNN imputation
dataset <- cached_block(
  file = "precalc_qmdiab.rds", 
  # checksum = "",
  block={
    
    result <- c()
    # load unprocessed QMdiab dataset, or check local file in current directory
    file <- "QMDiab_metabolomics_OrigScale.xlsx"
    load.web.file(
      url="https://ndownloader.figshare.com/files/10531339",
      md5sum = 'f61ad81cabaa3c3334b42ba1c46288a7',
      outfile = file
    )
    
    # loop over three fluids
    for (sheet in c("plasma","urine","saliva")) {
      
      # load metabolomics data
      raw <- read_excel(file, sheet = sheet)
      # split into assay, colData and rowData
      assay <- raw %>% dplyr::select(-`QMDiab-ID`, -AGE, -GENDER, -BMI, -ETHNICITY, -T2D)
      colData <- raw %>% dplyr::select(AGE, GENDER, BMI, ETHNICITY,T2D)
      rowData <- read_excel(file, sheet = sprintf("%s annotations",sheet), col_names = T) %>%
        dplyr::mutate(rn=BIOCHEMICAL, name=BIOCHEMICAL) %>%
        tibble::column_to_rownames("rn")
      
      # run preprocessing
      D <- mapletify_data(assay=assay,colData=colData,rowData=rowData) %>%
        maplet::mt_anno_mutate(anno_type = "features", col_name = "HMDB", term=HMDb_ID) %>%
        maplet::mt_anno_mutate(anno_type = "samples", col_name = "group", term=as.factor(T2D)) %>%
        fun_preprocessing()
      
      data <- list()
      
      # save preprocessed data
      data$X <- D$norm %>% assay %>% t %>% as.data.frame()
      # save minimum imputed data
      data$X.minimp <- D$minimp %>% assay %>% t %>% as.data.frame()
      # save knn imputed data
      data$X.knnimp <- D$knnimp %>% assay %>% t %>% as.data.frame()
      
      # save covariates
      data$Covars <- D$norm %>% 
        colData() %>% as.data.frame() %>% 
        dplyr::select(AGE, GENDER, BMI, ETHNICITY)
      # save group annotations
      data$groups <- D$norm %>% 
        colData() %>% as.data.frame() %>% 
        dplyr::pull(T2D)
      # save metabolite annotations
      data$anno <- D$norm %>% rowData() %>% as.data.frame() %>%
        dplyr::select(BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY, KEGG, HMDB)
      
      # store
      result[[paste0("qmdiab.", sheet)]] <- data
      
    }
    # return
    result
  })

# add to list of datasets
datasets <- c(datasets, dataset)



#### Breast Cancer Data (Terunuma et al., 2014) ----

# Original Publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3871244/

# cached block, to avoid recomputing of KNN imputing
dataset <- cached_block(
  file = "precalc_breast.rds", 
  # checksum = "",
  block={
    
    result <- c()
    # load unprocessed dataset, or check local file in current directory
    file <- "JCI71180sd2.xlsx"
    load.web.file(
      url="https://dm5migu4zj3pb.cloudfront.net/manuscripts/71000/71180/JCI71180sd2.xlsx",
      md5sum = 'd2c76018e6025479632e8a951d4fb4b0',
      outfile = file
    )
    
    sheet <- "OrigData"
    # load metabolomics data
    raw = read_excel(path=file, sheet=sheet, col_names = F)
    raw[1,1] <- NA
    
    # split into assay, colData and rowData
    assay <- raw[21:dim(raw)[1],11:dim(raw)[2]] %>% t %>% as.data.frame()
    assay %<>% apply(2, as.numeric)
    colnames(assay) <- raw[21:(20+ncol(assay)),1] %>% unlist() %>% as.vector
    rownames(assay) <- raw[2,11:(10+nrow(assay))] %>% unlist() %>% as.vector
    
    colData <- raw[4:20,11:dim(raw)[2]] %>% t() %>% as.data.frame()
    colnames(colData) <- raw[4:20,10] %>% unlist() %>% make.names
    rownames(colData) <- rownames(assay)
    colData %<>% 
      dplyr::rename(ETHNICITY=RACE.ETHNICITY) %>%
      dplyr::mutate(groups=ifelse(TISSUE.TYPE %in% c("POS TUMOR","NEG TUMOR"),"TUMOR","NORMAL") %>% as.factor)
    
    rowData <- raw[21:dim(raw)[1],1:10]
    colnames(rowData) <- raw[20,1:10]
    rowData %<>%
      dplyr::rename(HMDB=`HMDB_ID\\TISSUE WEIGHT mg`) %>%
      dplyr::rename(KEGG=KEGG_ID)  %>%
      dplyr::mutate(rn=BIOCHEMICAL, name=BIOCHEMICAL) %>%
      tibble::column_to_rownames("rn")
    
    # run preprocessing
    D <- mapletify_data(assay=assay,colData=colData,rowData=rowData) %>%
      fun_preprocessing()
    
    data <- list()
    
    # save preprocessed data
    data$X <- D$norm %>% assay %>% t %>% as.data.frame()
    # save minimum imputed data
    data$X.minimp <- D$minimp %>% assay %>% t %>% as.data.frame()
    # save knn imputed data
    data$X.knnimp <- D$knnimp %>% assay %>% t %>% as.data.frame()
    
    # save covariates
    data$Covars <- D$norm %>% 
      colData() %>% as.data.frame() %>% 
      dplyr::select(AGE, ETHNICITY)
    # save group annotations
    data$groups <- D$norm %>% 
      colData() %>% as.data.frame() %>% 
      dplyr::pull(groups)
    # save metabolite annotations
    data$anno <- D$norm %>% rowData() %>% as.data.frame() %>%
      dplyr::select(BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY, KEGG, HMDB)
    
    # store
    setNames(list(data), "brca")
    result[["brca"]] <- data
    
    # return
    result
    
  })

# add to list of datasets
datasets <- c(datasets, dataset)



#### Kidney Cancer Data (Hakimi et al., 2016) ----

# Original Publication: https://www.cell.com/cancer-cell/fulltext/S1535-6108(15)00468-7#secsectitle0145

# cached block, to avoid recomputing of KNN imputing
dataset <- cached_block(
  file = "precalc_kidney.rds", 
  # checksum = "",
  block={
    
    result <- c()
    # load unprocessed dataset, or check local file in current directory
    file <- "mmc2.xlsx"
    load.web.file(
      url="https://www.cell.com/cms/10.1016/j.ccell.2015.12.004/attachment/59785359-a5a4-4065-89ae-6a72abc88084/mmc2.xlsx",
      md5sum = '9408e74e2a92691433c16d7a850f82c5',
      outfile = file
    )
    
    sheet <- "Raw Data"
    # load metabolomics data
    raw = read_excel(path=file, sheet=sheet, col_names = F)
    
    # split into assay, colData and rowData
    assay <- raw[13:dim(raw)[1],13:dim(raw)[2]] %>% t() %>% as.data.frame()
    assay %<>% apply(2, as.numeric)
    colnames(assay) <- raw[13:(12+ncol(assay)),2] %>% unlist() %>% as.vector()
    rownames(assay) <- raw[12,13:(12+nrow(assay))] %>% unlist() %>% as.vector()
    
    colData <- raw[2:11,13:dim(raw)[2]] %>% t() %>% as.data.frame()
    colnames(colData) <- raw[2:11,12] %>% unlist() %>% make.names
    rownames(colData) <- rownames(assay)
    colData %<>%
      dplyr::rename(AGE=AGE.AT.SURGERY,
                    ETHNICITY=RACE) %>%
      dplyr::mutate(groups=ifelse(TISSUE.TYPE == "T", "TUMOR","NORMAL") %>% as.factor)
    
    rowData <- raw[13:dim(raw)[1],1:12]
    colnames(rowData) <- raw[12,1:12]
    rowData %<>% 
      dplyr::rename(BIOCHEMICAL=`BIOCHEMICAL NAME`) %>%
      dplyr::rename(HMDB=`Group HMDB`) %>%
      dplyr::rename(SUPER_PATHWAY=`SUPER PATHWAY`) %>%
      dplyr::rename(SUB_PATHWAY=`SUB PATHWAY`) %>%
      dplyr::mutate(rn=BIOCHEMICAL, name=BIOCHEMICAL) %>%
      tibble::column_to_rownames("rn")
    
    # run preprocessing
    D <- mapletify_data(assay=assay,colData=colData,rowData=rowData) %>%
      fun_preprocessing()
    
    data <- list()
    
    # save preprocessed data
    data$X <- D$norm %>% assay %>% t %>% as.data.frame()
    # save minimum imputed data
    data$X.minimp <- D$minimp %>% assay %>% t %>% as.data.frame()
    # save knn imputed data
    data$X.knnimp <- D$knnimp %>% assay %>% t %>% as.data.frame()
    
    # save covariates
    data$Covars <- D$norm %>% 
      colData() %>% as.data.frame() %>% 
      dplyr::select(AGE, GENDER, ETHNICITY)
    # save group annotations
    data$groups <- D$norm %>% 
      colData() %>% as.data.frame() %>% 
      dplyr::pull(groups)
    # save metabolite annotations
    data$anno <- D$norm %>% rowData() %>% as.data.frame() %>%
      dplyr::select(BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY, KEGG, HMDB)
    
    # return
    setNames(list(data), "rcc")
    
  })

# add to list of datasets
datasets <- c(datasets, dataset)



# HAPO dataset ----
# M Scholtens et al., https://doi.org/10.2337/dc13-0989

# data comes from metabomxtr
# package if not installed, these lines need to be executed
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("metabomxtr")

data(metabdata, package = "metabomxtr")

datasets$metabomxtr <- list(X = metabdata[,-seq(10)], Covars = NULL, groups = metabdata$PHENO, X.minimp = NULL, X.knnimp = NULL )
datasets$metabomxtr$X.minimp <- datasets$metabomxtr$X %>% apply(2,function(x){x[is.na(x)] = min(x,na.rm = T);x})
datasets$metabomxtr$X.knnimp <- datasets$metabomxtr$X %>% maplet:::imputeKNN(methods = "knn.obs.euc.sel", K = 10)

# save prepared data 
save(file = "metabolomics_datasets.Rdata", datasets)



# status
cat("\n\n.Rdata files saved to current working directory.\n\n")

