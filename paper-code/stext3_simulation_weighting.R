# This script generates the figure in Supplementary Text 3


rm(list = ls())
library(survival)
library(magrittr)
library(ggplot2)

n = 1e4
set.seed(42)
x = rnorm(n) # metabolite
z = rnorm(n) # noise
z = z[order(x)]
x = sort(x)

# pick the parameter such that desired concordance is obtained: 
#     d = concordance(y~x), where y = x + nfs*z
# noise factors
nfs = seq(sqrt(0.5),sqrt(6.5),length.out = 601)^2

# here I find corresponding noise factor for the sequence of ds 
#    ds = 0.55, 0.6, 0.65...
## find noise factors for equidistant ds 
# ds<-
# sapply(nfs, function(noise_factor){
#   y = z*noise_factor + x
#   concordance(x~y)$concordance
# })

# you can run the code above and confirm, these settings give 
# desired concordance values
ii = c( 573, 336, 224, 151, 95, 48, 3)
nfs = nfs[ii] 


# function to calculate optimal weight for x,y
#             given the missingness percentage
ff <- function(x, y, pcth){
  th = x[length(x)*pcth] - .Machine$double.eps
  yh = 1*(x > th)
  xh = max(x) -x + 0.1
  xh[yh==0] = max(xh)
  S = Surv(xh, yh)
  
  # concordance for pi_1
  d1 = concordance(x~y, data.frame(x,y)[yh==1,])$concordance
  # concordance for pi_2
  d2 = concordance(x~y, data.frame(x,y)[yh==0,])$concordance
  # concordance for bridge pairs
  db = concordance(yh~y)$concordance
  # concordance for left-censored formulation without weighting
  dS = concordance(S~y, reverse = T)$concordance
  # true concordance
  d = concordance(x~y)$concordance
  
  # optimal weight, formula derived in the paper, see the paper
  w = (d-d1)/(db-d) * (db-dS)/(dS-d1)
  
  # if we apply optimal weighting, d==dw
  dw = if(w>0) concordance(S~y, reverse = T, weights = c(w,1)[yh+1])$concordance else NA
  
  c(w=w, m = pcth, dw=dw, d = d, dS = dS, d1=d1, d2=d2, db=db)
}

# missingness percentage
pcs = seq(0.1, 0.9, by = 0.05)
pcs = structure(pcs, names = pcs)

re = lapply( pcs, function(pcth) 
  data.frame( t(sapply(nfs, function(k) ff(x, z*k + x, pcth))), miss = pcth )
) %>% do.call(what = rbind)


# 500 x 500
ggplot(re, aes(x = miss, y= w, color = factor(round(d,2)))) + 
  # geom_line( fill = NA ) +
  geom_smooth( aes(fill =  factor(round(d,2)) ), alpha = 0.25, lwd = 0.5) + 
  geom_point()+
  scale_color_manual(values = rev(gray.colors(7))) + 
  scale_fill_manual(values = rev(gray.colors(7))) + 
  theme_classic() +
  geom_abline(slope = -1, intercept = 1,  lty = 2, color = "red") +
  labs(y="weight", x = "missingness ratio(m)", color = "concordance", fill = "concordance") +
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  coord_fixed() +
  annotate(geom = "label", x = 0.725, y = 0.65, label = 'proposed weight = 1-m', size = 3.5) +
  geom_curve(aes(x = 0.724, y = 0.6, xend = 0.63, yend = 0.38), lty = 1,
             arrow = arrow(length = unit(0.025, "npc"),type = "closed"), curvature = -0.2, lwd = 0.5)

ggsave(filename = 'ST3_weighting_simulation.pdf', width = 5, height = 5)



# status
cat("\nFigure saved as PDF in current working directory\n")

