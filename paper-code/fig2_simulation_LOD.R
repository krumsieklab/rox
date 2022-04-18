# This script generates Figure 2
#
# The execution takes about 50s on an Apple M1 Pro


rm(list= ls())

library(magrittr)
library(survival)
library(ggplot2)
library(patchwork)

# rox 
library(rox)

# run simulation if results have not been saved before ----

# utility functions for simulation
source('utils_simulation.R')

# sample size for simulation 
n = 10000

# missingness percentages
ps = seq(0, 0.9, by = 0.025)

# random data
set.seed(42)
x = rnorm(n) # metabolite
z = rnorm(n) # noise

# pick noise_factor such that desired concordance is obtained: 
#     d = concordance(y~x), where y = x + nfs*z
# true concordances for which noise_factor to be found 
ds = seq(0.55, 0.85, by = 0.05)
nfs = sapply(ds, find_noise_factor)
 
# run simulation analysis
df<-
  lapply(nfs, function(noise_factor){
    # simulated outcome
    y = z*noise_factor + x
    # true concordance
    d = round(concordance(x~y)$concordance, 2)
    re <- sapply(ps, function(p){
      # induce missigness 
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
    }) %>% t 
    data.frame( ps = ps, d = d, noise_factor = noise_factor, re)
  }) %>% do.call(what = rbind)



# Figure 2 ----------------------------------------------------------------

# melt the data.frame for ggplot
mdf = reshape2::melt( df, id = c("ps", "d", "noise_factor"))

# ground-truth concordances
ds = mdf$d %>% unique %>% sort
mdf$variable %<>% factor

# name the names 
mdf$variable = factor(mdf$variable, 
                      levels = c("ds", "d1", "dw"), 
                      labels = c("min-imp", "CCA", "rox"))
# keep only some ds 
mdf = mdf %>% dplyr::filter(round(d,2) %in% c(0.55, 0.65, 0.75, 0.85))

color_codes = c("steelblue","orange", "#45A445", "red")

# ggplot
gg_rox =
  ggplot(mdf, aes(ps*100, value, color = variable))+
  geom_line(aes(y=d), color= "black", lty = 2)+
  geom_line()+
  facet_wrap(~paste0("d=", round(d,2)), nrow = 1) + 
  labs(color = "", x= "missingness %", y = "concordance") + 
  theme_minimal() + 
  theme(legend.position = "bottom", 
        strip.background = element_rect(fill = "wheat", color = "black"),
        panel.background = element_rect(colour = "black") ) +
  scale_color_manual(values = color_codes)

# LOD-based missingness plot
gg_data = 
  ggplot(data.frame(x = c(-4, 4)), aes(x, fill = x<0)) +
  stat_function(fun = dnorm, xlim = c(-4,0), geom = "area", 
                alpha = 0.7, aes(fill = "missing")) +
  stat_function(fun = dnorm, xlim = c(0,4), geom = "area", 
                alpha = 0.7, aes(fill = "measured")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "wheat", color = "black"),
        panel.background = element_rect(colour = "black") ) +
  labs(fill = "", y = "density") +
  facet_grid(.~("%50 missingness with strict-LOD"), scales = "free") +
  scale_fill_manual(values = c("black","gray")) + 
  guides(fill = guide_legend(reverse = T))


# generate, 900 X 350
gg_data + 
  (gg_rox + 
     aes(lty = "ground truth") + 
     labs(lty = NULL) +
     guides(
       linetype = guide_legend(override.aes = list(
         shape = c(NA),
         linetype = c("dashed")
       ), nrow = 1, byrow = TRUE)
     )
  )+ 
  patchwork::plot_layout(nrow = 1, widths = c(3,8)) +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = 'Fig2_simulation_LOD.pdf', width = 9, height = 3.5)

# status
cat("\nFigure saved as PDF in current working directory\n")



