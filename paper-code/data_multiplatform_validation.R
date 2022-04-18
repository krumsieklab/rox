# This script preprocesses and prepares the data for the multi-platform validation analysis
# It requires Supplementary Data 1 as "multi_platform_validation_supplement.xlsx" in the working directory
#
# The execution takes about 13min on an Apple M1 Pro

rm(list = ls())
library(SummarizedExperiment)
library(magrittr)
# maplet
if (!("maplet" %in% rownames(installed.packages()))) stop("This script needs the maplet package. devtools::install_github(repo='krumsieklab/maplet', subdir='maplet')")
library(maplet)


# metabolomics measured with HD2 platform
hd2 = readxl::read_xlsx("multi_platform_validation_supplement.xlsx", sheet = "HD2") %>% {
  a = as.matrix(.[,-1]);
  rownames(a) =unlist(.[,1]);
  a
}

# metabolomics measured with HD4 platform
hd4 = readxl::read_xlsx("multi_platform_validation_supplement.xlsx", sheet = "HD4") %>% {
  a = as.matrix(.[,-1]);
  rownames(a) =unlist(.[,1]);
  a
}

# phenotype data 
phenotypes = readxl::read_xlsx("multi_platform_validation_supplement.xlsx", 
                               sheet = "phenotypes") %>% as.data.frame()

# convert to SE objects
S2 = list( 
  hdf = SummarizedExperiment(assays = list(t(hd4)), colData = phenotypes[,-1]),
  pm = SummarizedExperiment(assays = list(t(hd2)), colData = phenotypes[,-1])
)

# quotient normalization and log scaling
library(maplet)
S2 = lapply(S2, function(D){
  
  # change 0s to NAs
  assay(D)[assay(D) == 0] <- NA
  D %>%
    mt_pre_filter_missingness(feat_max = 0.8) %>%
    mt_pre_norm_quot(feat_max = 0.2) %>%
    mt_pre_trans_log() 
})

# raw data
dt_raw = lapply(S2, function(x) t(assay(x)))
i = lapply(dt_raw, colnames) %>% {intersect(.[[1]], .[[2]] )}
j = lapply(dt_raw, rownames) %>% {intersect(.[[1]], .[[2]] )}
dt_raw = lapply(dt_raw, `[`,j,i)

# minimum imputation
fimp_min <- function(x) apply(x, 2, function(x){
  x[is.na(x)] = min(x, na.rm = T)
  x
})
dt_minp = lapply(dt_raw, function(x) fimp_min(x))

# knn imputation
fimp_knn <- maplet:::imputeKNN
dt_knn = lapply(dt_raw, function(x) fimp_knn(x))

# outcomes of interest 
df_y = as.data.frame( colData(S2$hdf)[j,] )
df_y$AGE = scale(df_y$AGE)
df_y$BMI = scale(df_y$BMI)


# find metabolites to be considered ---------------------
# here we decide between Fully Quantified (FQ) and Partially Missing (PM) metabolites

# summary of cross missingness pattern
sm = 
  seq(ncol(dt_raw$hdf)) %>% structure(., names = colnames(dt_raw$hdf)) %>% 
  sapply(function(i){
    x1 = dt_raw$hdf[,i]
    x2 = dt_raw$pm[,i]
    tb = table(c(T,F,!is.na(x1)), c(T,F,!is.na(x2)))
    diag(tb) = diag(tb)-1
    tb
    c(complete = tb[2,2],  
      miss_max = max(diag(tb[,2:1])), 
      miss_min = min(diag(tb[,2:1]))+tb[1,1] )
  }) %>% t %>% data.frame()

# cross missingness, effective missingness percentage
xmissing = sm[,2]/rowSums(sm)
names(xmissing) = colnames(dt_raw$hdf)

# all missing, missing combined 
amissing = (sm[,2]+sm[,3])/ rowSums(sm)

# lost anyway 
lostanyway = sm[,"miss_min"]/ rowSums(sm)

# metabolites to be considered
jinds = which( xmissing>=0 & lostanyway < 0.05 )
cut(xmissing[jinds], breaks = c(-1, 0, 0.05, 0.1, 0.15 ,0.2, 1)) %>% table 

# keep only those metabolites that are considered
dt_raw = lapply(dt_raw, `[`, ,jinds)
dt_minp = lapply(dt_minp, `[`, ,jinds)
dt_knn = lapply(dt_knn, `[`, ,jinds)

# save the data for subsequent analysis
save(file = "qmdiab_plasma_hd2_hd4.Rdata", dt_raw, dt_minp, dt_knn, df_y)

# status
cat("\n\n.Rdata files saved to current working directory.\n\n")
