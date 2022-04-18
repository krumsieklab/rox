# This script generates Supplementary Figure 1
#
# The execution takes about 7min on an Apple M1 Pro

rm(list = ls())
library(magrittr)
library(survival)
library(ggplot2)
library(rox)


# run simulation  -------------------
# utility functions for simulation
source('utils_simulation.R')

# sample size for population
n = 1000000

# generated random data
set.seed(42)
x0 = rnorm(n)   # metabolite
z0 = rnorm(n)   # noise 

# pick noise_factor such that desired concordance is obtained: 
#     d = concordance(y~x), where y = x + nfs*z
# true concordances for which noise_factor to be found 
ds = c(0.6, 0.7, 0.8)
nfs = sapply(ds, find_noise_factor)
             
# ground truth concordances
ds<-
  sapply(nfs, function(noise_factor){
    y0 = z0*noise_factor + x0
    d = concordance(x0~y0)$concordance
  })

# small sample size 
n1 = 100
# missigness percentages
ps = c(0.1, 0.3, 0.5, 0.7, 0.9)

# number of small sample size samples
K = 1000 # replicate


# LOD simulation ----------------------------------------------------------
  
set.seed(42)
df_lod<-
  lapply(seq(nfs), function(i){
    print(i)
    noise_factor = nfs[i]
    # ground truth 
    d = ds[i]
    
    re <- replicate(K, sapply(ps, function(p){
      # sample randomly from the population
      inds = sample(n,n1)
      x = x0[inds]
      y = z0[inds]*noise_factor + x
      # induce missingness
      xh = fnaify(x,p)
      inds = !is.na(xh)
      
      # CCA
      d1 = concordance(xh[inds]~y[inds])$concordance
      # rox
      dw = unname(rox(xh, y)$stats["d"])
      # min imputation
      xminp = xh %>% {.[is.na(.)]=min(.,na.rm = T);.}
      dm = concordance(xminp~y)$concordance
      
      c(d1=d1, dw=dw, ds=dm)
    }) )
    list( ps = ps, d = d, noise_factor = noise_factor, re = re)
  })
  
  
# pLOD simulation ---------------------------------------------------------

# pLOD values for simulating probabilistic LOD
# range LODness from strict-LOD to missingness-at-random(MARS)
plods = c( seq(1, 0.6, by = -0.05), 0.4, 0.2, 0 )

set.seed(42)
df_plod<-
  lapply(seq(nfs), function(i){
    print(i)
    noise_factor = nfs[i]
    # ground truth 
    d = ds[i]
    
    re <- replicate(K, sapply(plods, function(plod){
      # sample randomly from the population
      jinds = sample(n,n1)
      x = x0[jinds]
      y = z0[jinds]*noise_factor + x
      # induce missingness
      xh = fnaify_plod(x, plod)
      inds = !is.na(xh)
      
      # CCA
      d1 = concordance(xh[inds]~y[inds])$concordance
      # rox
      dw = unname(rox(xh, y)$stats["d"])
      # min imputation
      xminp = xh %>% {.[is.na(.)]=min(.,na.rm = T);.}
      dm = concordance(xminp~y)$concordance
      
      c(d1=d1, dw=dw, ds=dm)
    }))
    list(ss =plods, d = d, noise_factor = noise_factor, re = re)
  })




# generate the figure ----------------------------------------------------

# compiles the data into ggplot frendly format
get_data <- function(df_){
  mm =
    lapply(df_, function(x){
      mm = x$re[c("d1","dw","ds"),,] %>% reshape2::melt()
      mm$Var3 = NULL
      colnames(mm) = c("model", "miss", "d")
      mm$d0 = x$d
      mm
    }) %>% do.call(what = rbind)
  
  mm$d0= round(mm$d0, digits = 2)
  levels(mm$model) = c(d1 = "CCA", dw = "rox", ds = "min-imp")[levels(mm$model)]
  if(any(mm$miss>1)) mm$miss = df_[[1]][[1]][mm$miss]
  mm
}

# plot LOD results
mm = get_data(df_lod)
sm = mm %>% dplyr::group_by(d0, miss, model) %>% dplyr::summarise(d = median(d))

gg1 = 
  ggplot(mm, aes(y = d,color = model,x =  factor(miss*100), fill = model)) +
  geom_hline(aes(yintercept = d0), color = "gray30", lty = 2) +
  geom_path(data = sm, aes(group = model)) +
  geom_boxplot(width = 0.2, outlier.colour = NA, position=position_dodge(0.4), alpha = 0.25) +
  facet_grid(.~factor(d0, levels = c(0.6,0.7, 0.8), 
                      labels = c("d=0.6", "d=0.7", "d=0.8"))) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "wheat"), 
        panel.background = element_rect()) +
  labs(y= "concordance", x = "missingness %", color = "", fill = "") + 
  coord_fixed(ratio = 4.5) 


# plot prob-LOD results
mm = get_data(df_plod)
mm = mm %>% dplyr::filter(miss %in% ((0:5)/5))
sm = mm %>% dplyr::group_by(d0, miss, model) %>% dplyr::summarise(d = median(d))

gg2=
  ggplot(mm, aes(y = d,color = model,x =  factor(miss), fill = model)) +
  geom_hline(aes(yintercept = d0), color = "gray30", lty = 2) +
  geom_path(data = sm, aes(group = model)) +
  geom_boxplot(width = 0.2, outlier.colour = NA, position=position_dodge(0.4), alpha = 0.25) +
  facet_grid(.~factor(d0, levels = c(0.6,0.7, 0.8), 
                      labels = c("d=0.6", "d=0.7", "d=0.8"))) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "wheat"), 
        panel.background = element_rect()) +
  labs(y= "concordance", x = "pLOD", color = "", fill = "") + 
  coord_fixed(ratio = 4.5 * 2) 


color_codes = c("orange", "#45A445","steelblue", "red")

gg1 = gg1 + 
  scale_color_manual(values = color_codes) +
  scale_fill_manual(values = color_codes)


gg2 = gg2 + 
  scale_color_manual(values = color_codes) +
  scale_fill_manual(values = color_codes)

# 900 x600
gg1/gg2

ggsave(filename = 'SFig1_simulation_small_sample_size.pdf', width = 9, height = 6)



# status
cat("\nFigure saved as PDF in current working directory\n")





