# This script generates Supplementary Figure 3

rm(list = ls())
library(magrittr)
library(survival)
library(ggplot2)
library(rox)

# keep only raw data with missing values
load('qmdiab_plasma_hd2_hd4.Rdata')
rm(list = setdiff(ls(), "dt_raw"))

# calculate effective sample size for each dataset for each metabolites [@mubu, not clear what "effective sample size" means]
# the data having more samples measures will be accepted as scores 
# the other data will be used to define binary classes
a1=
  seq(ncol(dt_raw$hdf)) %>% 
  sapply(function(i) 
    concordance(as.numeric(is.na(dt_raw$hdf[,i]))~dt_raw$pm[,i])$count[c(1,2)] %>% 
      sum %>% {sqrt(2*.)}
  )

a2=
  seq(ncol(dt_raw$hdf)) %>% 
  sapply(function(i) 
    concordance(as.numeric(is.na(dt_raw$pm[,i]))~dt_raw$hdf[,i])$count[c(1,2)] %>% 
      sum %>% {sqrt(2*.)}
  )


# how many data points ultimately each test will be based on [@mubu, fix english, not a sentence]
# which platform will be used for class definition missing yes, no, 
# and which platfrom will be used as scores in AUC calculation
apply( cbind(a1, a2), 1, max )
j = apply( cbind(a1, a2), 1, which.max )

# calculate AUCs 
ds =
  seq(ncol(dt_raw$hdf)) %>% 
  sapply(function(i){
    cfit = concordance(as.numeric(is.na(dt_raw[[j[i]]][,i]))~dt_raw[[3-j[i]]][,i])
    cfit$concordance
  })


# simple version 
data.frame(ds=1-ds) %>% na.omit %>% 
  ggplot(aes(y=ds, x = "pLODh")) +
  geom_hline(yintercept = 0.5, color = "red", lty = 3) +
  geom_violin(fill = NA,draw_quantiles = c(0.25, 0.5, 0.75), width = 1.1, lwd = 0.5) +
  ggbeeswarm::geom_quasirandom(color  = "blue") +
  #geom_jitter(aes(color = z), width = 0.25, pch = 21) + 
  theme_minimal() +
  scale_color_manual(values = c("blue", "gray", "brown")) +
  labs(x = "", color = "LODness",
       #y = latex2exp::TeX(r"($P(x_i < x_j | x_i = NA, x_j \neq NA)$)")) +
       y = latex2exp::TeX(r"($AUC_{LOD}$)")) +
  theme(axis.text.x = element_blank()) + 
  ylim(c(0,1))


# 500 x 600
ggsave(filename = 'SFig3_cross_lod_estimation.pdf', width = 5, height = 6)



# status
cat("\nFigure saved as PDF in current working directory\n")






