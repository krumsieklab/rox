# This script generates Figure 3

rm(list= ls())
library(magrittr)
library(survival)
library(ggplot2)
library(patchwork)

# rox 
library(rox)

# run simulation ----

# utility functions for simulation
source('utils_simulation.R')

# sample size for simulation 
n = 10000

# example simulations for continuous outcome
set.seed(42)
x = rnorm(n) # metabolite
z = rnorm(n) # noise

# pick noise_factor such that desired concordance is obtained: 
#     d = concordance(y~x), where y = x + nfs*z
# true concordances for which noise_factor to be found 
ds = seq(0.55, 0.85, by = 0.05)
nfs = sapply(ds, find_noise_factor)

# desired pLOD for simulating probabilistic LOD
# range the LODness from strict-LOD to missingness-at-random(MARS)
plods = c( seq(1, 0.6, by = -0.05), 0.4, 0.2, 0 )

df<-
  lapply(nfs, function(noise_factor){
    # simulated outcome
    y = z*noise_factor + x
    # true concordance
    d = round(concordance(x~y)$concordance,2)
    
    re <- sapply(plods, function(plod){
      # induce missingness with probabilistic LOD
      xh = fnaify_plod(x, plod)
      inds = !is.na(xh)
      
      #CCA
      d1 = concordance(xh[inds]~y[inds])$concordance
      # rox
      dw = unname(rox(xh, y)$stats["d"])
      # min imputation
      xminp = xh %>% {.[is.na(.)]=min(.,na.rm = T);.}
      dm = concordance(xminp~y)$concordance
      
      c(d1=d1, dw=dw, ds=dm)
    }) %>% t 
    data.frame(ss = plods, d = d, noise_factor = noise_factor, re)
    
  }) %>% do.call(what = rbind)

# data to exemplify probabilistic-LOD  
dxf<- 
  lapply(seq(1, 0, length.out = 6), function(plod){
    xh = fnaify_plod(x, plod)
    y = scale(z*nfs[3] + x)
    data.frame(x=x, y=y, xna = is.na(xh), s = as.character(plod) )
  }) %>% do.call(what = rbind)




# Figure 3 ----------------------------------------------------------------

# when dw differs from d1 
df$lwx = df$dw<df$d1

# melt the data for ggplot
mdf = reshape2::melt( df, id = c("ss", "d", "noise_factor", "lwx"))

# ground-truth concordances
ds = mdf$d %>% unique %>% sort 

mdf$variable %<>% factor
# line width
mdf$lsize = with(mdf, (!(lwx & variable == "d1") + (!lwx & variable == "dw"))>0)

# name the names 
mdf$variable = factor(mdf$variable, 
                      levels = c("ds", "d1", "dw"), 
                      labels = c("min-imp", "CCA", "rox"))
# keep only some d  
mdf = mdf %>% dplyr::filter(round(d,2) %in% c(0.55, 0.65, 0.75, 0.85))

color_codes = c("steelblue","orange", "#45A445", "red")

# generate plot
gg_rox = 
  ggplot(mdf ,aes(ss, value, color = variable))+
  geom_rect(data = df %>% dplyr::filter(dw==d1) %>% 
              dplyr::filter(round(d,2) %in% c(0.55, 0.65, 0.75, 0.85)) %>%
              dplyr::group_by(d) %>% 
              dplyr::slice(which.max(ss)),
            aes(xmin=ss, xmax=1, ymin=-Inf, ymax=Inf), fill = "gold", inherit.aes = F, alpha = 0.17, color = "white") + 
  geom_line(aes(size = lsize))+
  geom_line(aes(y=d), color= "black", lty = 2, alpha = 0.5)+
  facet_grid(~paste0("d=", round(d,2))) + 
  labs(color = "", x= "pLOD", y = "concordance", size = "") + 
  theme_minimal() + 
  theme(legend.position = "right", 
        strip.background = element_rect(fill = "wheat", color = "black"),
        panel.background = element_rect(colour = "black") ) + 
  scale_size_manual(values = c(1.5,0.5)) +
  guides(size = "none")+
  theme(axis.text.x = element_text(size = 7)) +
  scale_x_continuous(breaks = (0:5)/5) +
  scale_color_manual(values = color_codes)

# pLOD figure 
give_me_plod <- function(){
  # unifotm distribution
  # for [-4,4], perfect LOD to MARS
  # 1, 0.8, 0.6, 0.4, 0.2, 0
  s = c(1e-10, 0.58, 1.22, 2.22, 4.87, Inf)
  
  x = runif(1000,-4,4)
  df = sapply(s, function(i) plogis(-x,scale = i))
  colnames(df) <- paste0("s=", round(s,2))
  df = reshape2::melt(df)
  df$Var1 = rep(x, length(s))
  
  colnames(df)[1:2] <- c("x","s")
  
  df$s = factor(df$s, labels = c("1.0","0.8","0.6", "0.4","0.2","0.0"))
  
  ggplot(df, aes(x=x, y = value)) + 
    geom_line(aes(group = s, color = s), lwd = 5, alpha = 0.05) + 
    geom_line(aes(group = s, color = s)) + 
    theme_light() + 
    labs(y = "P(missing|x)",
         color = "pLOD") + 
    # ggtitle("missing values are generated probabilistically based on P(LOD)") + 
    scale_color_viridis_d(end = 0.9, option = "E") +
    theme(axis.text = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8))
  
}

# generate data distributions
gg_data = 
  ggplot(dxf %>% dplyr::filter(! s %in% c("0", "1")), 
         aes(x,fill = factor(xna, labels = c("measured", "missing")) )) +
  geom_density(aes(y=..density..),color = NA, alpha = 0.7, trim = T, bw = 0.25) +
  facet_wrap(.~paste0("pLOD = ",s), ncol = 2) +
  scale_fill_manual(values = c("black","gray")) + #
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "wheat", color = "black"),
        panel.background = element_rect(colour = "black") ) +
  labs(fill = "")


# assemble and export figure
gg1 = give_me_plod()

gg2 = gg_data +
  theme(legend.position = "right") + 
  guides(fill = guide_legend(reverse = T))

gg3 = gg_rox + 
  aes(lty = "ground truth") + 
  labs(lty = NA) + 
  theme(legend.title = element_blank()) +
  guides(
    linetype = guide_legend(override.aes = list(
      shape = c(NA),
      linetype = c("dashed")
    ), nrow = 1, byrow = TRUE)
  )

( gg1 + gg2 + theme(legend.position = "right") +
    patchwork::plot_layout(ncol = 2,nrow = 1, widths = c(1,2)) )/ ( gg3 ) + 
  patchwork::plot_layout(1,nrow = 2, heights = c(1,1.5)) +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = 'Fig3_simulation_probabilistic_LOD.pdf', width = 9, height = 6)

# status
cat("\nFigure saved as PDF in current working directory\n")

