# This script generates Supplementary Figure 2
#
# The execution takes about 4min on an Apple M1 Pro

rm(list = ls())
library(magrittr)
library(survival)
library(ggplot2)
library(rox)


# run simulation ----

# utility functions for simulation
source('utils_simulation.R')

# sample size for simulation
n = 10000

set.seed(42)
y = rnorm(n) # outcome
z = y + rnorm(n) # covariate
x =  y + z + 3*rnorm(n) # metabolite

# rox model
b0 = rox_mv(x,data.frame(z,y))
y = y*b0$stats_model['y', "B"]

# lm model
r0 = coef(lm(x~z+y))

# scale accordingly to have same beta
x = x/r0["y"]
lfit0 = lm(x~z+y)
dfit0 = rox_mv(x,data.frame(z,y))

# ground truth betas
rbind(
  lm = coef(lfit0)[-1],
  rox = dfit0$stats_model[, "B"]
) %>% round(2)


#  missingness percentages
mms = (1:9) / 10
# run the simulation
re_lod<-
  lapply(mms, function(ms){
    xi = x
    yi = y
    zi = z

    # induce missing values
    xh <- xminp <- fnaify(xi, ms, rand = F)

    # CCA
    inds = !is.na(xh)
    lfit_cca = lm(xi~zi+yi, data.frame(xi = xi[inds], zi = zi[inds], yi = yi[inds]))

    # min imputed lm fit
    xminp[is.na(xminp)] = min(xminp, na.rm = T)
    lfit_minimp = lm(xminp~zi+yi)

    # rox
    dfit_weighted =  rox_mv(xh, data.frame(z=zi, y=yi),)

    list(lfit_cca = lfit_cca, lfit_minimp = lfit_minimp, rox = dfit_weighted)
  })
names(re_lod) = mms

# save the results along with the data
nodes = list(data = data.frame(x,y,z),
             lfit0 = lfit0,
             dfit0 = dfit0,
             lod = re_lod)


# -----------------------------------------------------------------------
# probabilistic LOD results

set.seed(42)
# LODness
plods = (0:10) / 10
re_plod<-
  lapply(plods, function(pL){

    print(pL)
    xi = x
    yi = y
    zi = z

    # induce plod-based missingness
    xh <- xminp <- fnaify_plod(xi, pL)

    # CCA
    inds = !is.na(xh)
    lfit_cca = lm(xi~zi+yi, data.frame(xi = xi[inds], zi = zi[inds], yi = yi[inds]))

    # min imputed lm fit
    xminp[is.na(xminp)] = min(xminp, na.rm = T)
    lfit_minimp = lm(xminp~zi+yi)

    # weighted fits
    dfit_rox =  rox_mv(xh, data.frame(z=zi, y=yi))

    list(lfit_cca = lfit_cca, lfit_minimp = lfit_minimp, rox = dfit_rox)
  })
names(re_plod) = plods
nodes$plod = re_plod



# Supplementary Figure 2 --------------------------------------------------

# ground truth
gt = rbind( nodes$lfit0$coefficients[-1],nodes$dfit0$stats_model[,"B"]) %>%
  round(2) %>% apply(2, unique)

gt = data.frame(value = gt, variable = names(gt))
gt = rbind(gt, data.frame(value = nodes$dfit0$stats_overall["d"], variable = "d"))

# to obtain beta coefficients
f <- function(fit){
  if(class(fit)[1] == "roxfit")
    return( c(fit$stats_model[,"B"],fit$stats_overall["d"]) )
  c(fit$coefficients[-1], concordance(fit)$concordance)
}

# collect strict LOD results
df<-
  names(nodes$lod) %>% lapply(function(i){
    re = sapply(nodes$lod[[i]],f)
    data.frame(t(re), model = colnames(re), miss = as.numeric(i))
  }) %>% do.call(what = rbind)

colnames(df)[1:3] = c("z", "y", "d")

mdf = reshape2::melt(df, id = c("model", "miss"))
mdf$model = factor(mdf$model, levels = c("lfit_cca", "lfit_minimp", "rox"),
                   labels = c("lm + CCA", "lm + min-imp", "rox") )

color_codes = c("orange", "steelblue", "#45A445", "red")

# beta coefficients
gg_b=
  ggplot(mdf %>% dplyr::filter(variable != "d"), aes(x=miss*100, y = value, color = model)) +
  geom_hline(data = gt %>% dplyr::filter(variable != "d"), lty = 2,
             aes(yintercept = value), color = "black") +
  geom_line() +
  facet_grid(variable~., scales = "free_y") +
  geom_point(size = 2.5) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  labs(y = latex2exp::TeX("$\\beta$"), x = "missigness %") +
  scale_color_manual(values = color_codes)

# overall model fit
gg_d =
  ggplot(mdf %>% dplyr::filter(variable == "d"), aes(x=miss*100, y = value, color = model)) +
  geom_hline(data = gt %>% dplyr::filter(variable == "d"), lty = 2,
             aes(yintercept = value), color = "black") +
  geom_line() +
  geom_point(size = 2.5) +
  facet_grid(("overall model fit")~("strict-LOD")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  labs(y = "concordance", x = "missigness %") +
  scale_color_manual(values = color_codes)



# collect probabilistic LOD results
df<-
  names(nodes$plod) %>% lapply(function(i){
    re = sapply(nodes$plod[[i]],f)
    data.frame(t(re), model = colnames(re), plod = as.numeric(i))
  }) %>% do.call(what = rbind)

colnames(df)[1:3] = c("z", "y", "d")

mdf = reshape2::melt(df, id = c("model", "plod"))
mdf$model = factor(mdf$model, levels = c("lfit_cca", "lfit_minimp", "rox"),
                   labels = c("lm + CCA", "lm + min-imp", "rox") )

# beta coefficients
gg_plod_b=
  ggplot(mdf %>% dplyr::filter(variable != "d"), aes(x=plod, y = value, color = model)) +
  geom_hline(data = gt %>% dplyr::filter(variable != "d"), lty = 2,
             aes(yintercept = value), color = "black") +
  geom_line() +
  facet_grid(variable~., scales = "free_y") +
  geom_point(size = 2.5) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  labs(y = latex2exp::TeX("$\\beta$"), x = "pLOD") +
  scale_color_manual(values = color_codes)

# overall model fit
gg_plod_d =
  ggplot(mdf %>% dplyr::filter(variable == "d"), aes(x=plod, y = value, color = model)) +
  geom_hline(data = gt %>% dplyr::filter(variable == "d"), lty = 2,
             aes(yintercept = value), color = "black") +
  geom_line() +
  geom_point(size = 2.5) +
  facet_grid(("overall model fit")~("probabilistic-LOD")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  labs(y = "concordance", x = "pLOD")+
  scale_color_manual(values = color_codes)

patch1 = patchwork::wrap_plots(
  gg_d + theme(
    legend.position = "n",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ),
  gg_b + labs(color = "") +
    theme(legend.position = "n"),
  heights = c(1, 2)
)

patch2 = patchwork::wrap_plots(
  gg_plod_d + theme(
    legend.position = "n",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ),
  gg_plod_b + labs(color = ""),
  heights = c(1, 2)
)

# export figure
patchwork::wrap_plots(patch1, patch2,widths = c(1,1), nrow = 1 )

ggsave(filename = 'SFig2_simulation_mv.pdf', width = 10, height = 7.5)


# status
cat("\nFigure saved as PDF in current working directory\n")





