# This script generates Figure 5 and supplementary Figure 4 and 5
#
# It uses "qmdiab_plasma_hd2_hd4.Rdata", which is generated by 
# `data_multiplatform_validation.R`
#
# The execution takes about 40s on an Apple M1 Pro



rm(list = ls())
library(rox)
library(magrittr)
library(ggplot2)

# load the data
load("qmdiab_plasma_hd2_hd4.Rdata")

# run multi-platform validation analysis 

# univariate rox results ------------------------------------------------
nodes = lapply( structure(seq(ncol(dt_raw[[1]])), names = colnames(dt_raw[[1]])), function(i){
  print(i)
  x1 = dt_raw$hdf[,i]
  x2 = dt_raw$pm[,i]
  m1 = is.na(x1)
  m2 = is.na(x2)
  
  # decide fully quantified (FQ), and partially missing (PM) platforms
  j = seq(2)
  if(sum(m1) > sum(m2)){
    m = m2
  }else{
    m = m1
    j = rev(j)
  }
  
  re<-
    list(
      rox =
        lapply(df_y, function(y)
          structure( list( rox(x1[!m], y[!m]), 
                           rox(x2[!m], y[!m]))[j], names = c('PM', 'FQ') )
        ),
      minp=
        lapply(df_y, function(y)
          structure( list( rox(dt_minp$hdf[!m,i], y[!m]), 
                           rox(dt_minp$pm[!m,i], y[!m]))[j], names = c('PM', 'FQ') )
        ),
      knn =
        lapply(df_y, function(y)
          structure( list( rox(dt_knn$hdf[!m,i], y[!m]), 
                           rox(dt_knn$pm[!m,i], y[!m]))[j], names = c('PM', 'FQ') )
        ),
      cca = 
        lapply(df_y, function(y)
          structure( list( rox(x1[!(m1|m)], y[!(m1|m)]), 
                           rox(x2[!(m2|m)], y[!(m2|m)]))[j], names = c('PM', 'FQ') )
        )
    )
  
  re$miss = max(sum(m1[!m]), sum(m2[!m]))/length(m[!m])
  re$fully_quantified = c("hd2",'hd4')[j[1]]
  return(re)
})

# which platform is fully quantified (reference)
f = sapply(nodes, `[[`, 'fully_quantified') 
# effective missingness
m = sapply(nodes,`[[`,'miss')
  
# collect the results
methods = c("rox", "minp", "knn", "cca") %>% {structure(., names = .)}
outcomes= c('SEX', 'AGE', 'BMI', 'DIAB') %>% {structure(., names = .)}
re = 
  lapply(outcomes, function(o) lapply(methods, function(w){
    x = sapply( nodes, function(nod) sapply( nod[[w]][[o]], function(x) unname(x$stats['d']))) %>%t
    data.frame(miss = m, x, outcome = o, method = w, reference = f)
  }))
  
# transforms results to data frame
get_df <- function(re, outcomes = c('SEX', 'AGE', 'BMI', 'DIAB'), 
                   methods = c("rox", "minp", "knn", "cca")){
  lapply(outcomes, function(i){
    ths = seq(0, 0.2, by = 0.05)
    data.frame(
      th = ths, 
      outcome = i, 
      sapply(methods, function(j){
        mm = re[[i]][[j]]
        sapply(ths, function(th) cor( mm[mm$miss>th,]$FQ, mm[mm$miss>th,]$PM) )
      })
    )
  })
}
  

# results: multivariable model  -----------------------------------------
nodes_mv = lapply( structure(seq(ncol(dt_raw[[1]])), names = colnames(dt_raw[[1]])), function(i){
  print(i)
  x1 = dt_raw$hdf[,i]
  x2 = dt_raw$pm[,i]
  m1 = is.na(x1)
  m2 = is.na(x2)
  
  # decide fully quantified (FQ), and partially missing (PM) platforms
  j = seq(2)
  if(sum(m1) > sum(m2)){
    m = m2
  }else{
    m = m1
    j = rev(j)
  }
  
  # function to fit linear models 
  fit_lm <- function(x, yy) lm(x~.+0, data.frame(x=scale(x), yy))
  
  re<-
    list(
      rox = structure( list( rox_mv(x1[!m], df_y[!m,]), 
                             rox_mv(x2[!m], df_y[!m,]))[j], names = c('PM', 'FQ') ),
      
      minp= structure( list( fit_lm(dt_minp$hdf[!m,i], df_y[!m,]), 
                             fit_lm(dt_minp$pm[!m,i], df_y[!m,]))[j], names = c('PM', 'FQ') ),
      
      knn = structure( list( fit_lm(dt_knn$hdf[!m,i], df_y[!m,]), 
                             fit_lm(dt_knn$pm[!m,i], df_y[!m,]))[j], names = c('PM', 'FQ') ),
      
      cca = structure( list( fit_lm(x1[!(m1|m)], df_y[!(m1|m),]), 
                             fit_lm(x2[!(m2|m)], df_y[!(m2|m),]))[j], names = c('PM', 'FQ') )
    )
  
  
  re$miss = max(sum(m1[!m]), sum(m2[!m]))/length(m[!m])
  re$fully_quantified = c("hd2",'hd4')[j[1]]
  return(re)
})
# get the coef of rox fit
coef.roxfit <- function(fit) fit$stats_model[,'B']

# in which platform fully quantified (reference)
f = sapply(nodes_mv, `[[`, 'fully_quantified') 
# effective missingness
m = sapply(nodes_mv,`[[`,'miss')

# collect the results
methods = c("rox", "minp", "knn", "cca") %>% {structure(., names = .)}
outcomes= c('SEX', 'AGE', 'BMI', 'DIAB') %>% {structure(., names = .)}
re_mv = 
  lapply(outcomes, function(o) lapply(methods, function(w){
    x = sapply( nodes_mv, function(nod) sapply( nod[[w]], function(x) unname(coef(x)[o]))) %>%t
    data.frame(miss = m, x, outcome = o, method = w, reference = f)
  }))



# generate figure 5 ----------------------------------------------------------------

color_codes = c(
  minp = "steelblue", 
  cca = "orange", 
  rox = "#45A445", #scales::muted("green",c = 70, l =60), 
  knn = "red"
)

# collect the results for plotting
df = get_df(re) %>% 
  do.call( what = rbind ) %>% 
  reshape2::melt(id  = c("outcome","th")) 
colnames(df) =  c('y', 'miss', 'model','value')
df$model = factor(df$model, levels = rev(levels(df$model)))

# results
gg0 = 
ggplot(df, aes(miss*100,value, color = model, fill = model)) +
  geom_line(alpha = 1) +
  geom_point(size = 2)+
  geom_point(data = df %>% dplyr::filter(model == "rox" & y == "AGE" & miss == 0.2),
             pch = 13, fill = "black", size = 4.7, color = "black" )+
  facet_grid(.~factor(y,levels = c("AGE", "SEX", "BMI", "DIAB"))) +
  labs(x = "missingness >%", y = "consistency") +
  theme_minimal(12) +
  theme(strip.text = element_text(size = 14),
        strip.text.y = element_blank(),
        strip.background.x = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  coord_fixed(ratio = 65) + 
  scale_color_manual(values = color_codes[levels(df$model)])


# an example scatter plot for rox
adf = re$AGE$rox
adf = adf[adf$miss>0.2,]
colnames(adf)[5] = "model"
adf =  cbind(adf, met = rownames(adf))
adf$met[rowSums( abs(adf[,c('FQ', 'PM')] - 0.5)>0.05 )<1] = NA

library(ggrepel)
gg_ex = 
ggplot(adf, aes(FQ, PM, color = model, fill = model )) +
  geom_vline(xintercept = 0.5, color = "gray", lty = 2) +
  geom_hline(yintercept = 0.5, color = "gray", lty = 2) +
  geom_abline(slope = 1, intercept = 0, color = "gray", lty = 1, size = 2, alpha = 0.5) +
  geom_smooth( alpha = 0.15, aes(fill =  model), method = lm)+ 
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  guides(color = "none", fill = "none") +
  theme_minimal(11) +
  theme(strip.text = element_text(size = 9),
        strip.background.x = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        panel.grid = element_blank(), axis.ticks = element_line() ) +
  # scale_x_continuous(breaks = c(0.4,0.5,0.6),limits = c(0.4,0.62)) +
  # scale_y_continuous(breaks = c(0.4,0.5,0.6),limits = c(0.4,0.62))+
  scale_fill_manual(values = color_codes['rox']) +
  scale_color_manual(values = color_codes['rox']) +
  ggrepel::geom_text_repel(aes(label = met), color = "black", size = 3, 
                           box.padding = 0.75, min.segment.length = 0)  + 
  coord_fixed(ratio = 1) 

cowplot::plot_grid(
  gg_ex, gg0 + labs(color = "", fill = ""), ncol = 2, rel_widths = c(1,3), 
  labels = c("A", "B") )
# 1200 * 400 
ggsave(filename = 'Fig5_multiplatform_validation.pdf', width = 13.5, height = 4.5)



# generate figure 5 related supplement -----------------------------------------------------

# collect to results for detailed scatter plot
adf = do.call( unlist(re, recursive = F), what = rbind ) %>% {.[.$miss> 0.2,]}
colnames(adf)[5] = "model"
adf$outcome = factor(adf$outcome, levels = c('AGE','SEX','BMI','DIAB'))
  
ggs = 
ggplot(adf, aes(FQ, PM, color = model, fill = model)) +
  geom_vline(xintercept = 0.5, color = "gray", lty = 2) +
  geom_hline(yintercept = 0.5, color = "gray", lty = 2) +
  geom_abline(slope = 1, intercept = 0, color = "gray", lty = 1, size = 2, alpha = 0.5) +
  geom_smooth( alpha = 0.15, aes(fill =  model), method = MASS::rlm)+ # or "gam"
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  # facet_wrap(~paste(y,model,sep = "-"),scales = "free", nrow = 4) + 
  facet_grid(outcome~model) + 
  labs(title = "univariate analysis", subtitle = "missingness > %20") +
  theme_minimal(11) +
  theme(strip.text = element_text(size = 9),
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        panel.grid = element_blank(), 
        axis.ticks = element_line() )+
  scale_color_manual(values = color_codes) +
  scale_fill_manual( values = color_codes) +
  theme(legend.position = 'n')
  
# add correlations to each panel
ggs = 
  ggs + geom_text( 
    data = gg0$data[gg0$data$miss==0.2,] %>% 
      {colnames(.)[1] = 'outcome'; .$outcome = factor(.$outcome, levels = c('AGE','SEX','BMI','DIAB'));.}, 
    aes(x = 0.57, y=0.41, label = paste0("cor = ",round(value,2))), size = 3, color = "black") + 
  coord_fixed(1)

# accuracy
ggs0 = 
  ggplot(df, aes(miss*100,value, color = model, fill = model)) +
  geom_line(alpha = 1) +
  geom_point(size = 2)+
  geom_point(data = df %>% dplyr::filter( miss == 0.2),
             pch = 1, fill = "black", size = 4, color = "black" )+
  facet_grid(.~factor(y,levels = c("AGE", "SEX", "BMI", "DIAB"))) +
  labs(x = "missingness > %", y = "consistency") +
  theme_minimal(11) +
  theme(strip.text = element_text(size = 14),
        strip.text.y = element_blank(),
        strip.background.x = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  coord_fixed( ratio = 65 ) + 
  scale_color_manual(values = color_codes[levels(df$model)])


library(patchwork)
(ggs0 + labs(color = "", fill = ""))/ggs + 
  patchwork::plot_layout(heights = c(1,2.5)) + 
  patchwork::plot_annotation(tag_levels = 'A')

# 13 * 10 
ggsave(filename = 'SFig4_multiplatform_validation_detailed.pdf', width = 9, height = 12.5)



# supplementary figure: multivariable model -------------------------------
# collect the data
df = get_df(re_mv) %>% 
  do.call( what = rbind ) %>% 
  reshape2::melt(id  = c("outcome","th")) 
colnames(df) =  c('y', 'miss', 'model','value')
df$model = factor(df$model, levels = rev(levels(df$model)))

# plot the results
gg1 = 
  ggplot(df, aes(miss*100,value, color = model, fill = model)) +
  geom_line(alpha = 1) +
  geom_point(size = 1.5)+
  geom_point(data = df %>% dplyr::filter(y == "BMI" & miss != 0),
             pch = 1, fill = "black", size = 3, color = "black" )+
  facet_grid(.~factor(y,levels = c("AGE", "SEX", "BMI", "DIAB"))) +
  labs(x = "missingness >%", y = "consistency") +
  theme_minimal(12) +
  theme(strip.text = element_text(size = 14),
        strip.text.y = element_blank(),
        strip.background.x = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        axis.ticks = element_line()) +
  coord_fixed(ratio = 65) + 
  scale_color_manual(values = color_codes[levels(df$model)])


# collect to results for detailed scatter plot
adf=
lapply(c(0.05, 0.1,0.15,0.2), function(th){
  adf = do.call( unlist(re_mv, recursive = F), what = rbind ) %>% {.[.$outcome=='BMI',]}
  colnames(adf)[5] = "model"
  adf = adf[adf$miss>th,]
  adf$th = th
  adf
} ) %>% do.call(what = rbind)
adf$th = factor(paste0('missingness > ', adf$th))

gg1s = 
ggplot(adf, aes(FQ, PM, color = model, fill = model)) +
  geom_abline(slope = 1, intercept = 0, color = "gray", lty = 1, size = 2, alpha = 0.5) +
  geom_point( shape = 1, alpha = 0.7) +
  geom_smooth( alpha = 0.15, aes(fill =  model), method =lm)+ # or "gam"
  facet_grid(model~th) + 
  labs(title = "BMI") +
  theme_minimal(11) +
  theme(strip.text = element_text(size = 9),
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(),
        panel.grid = element_blank(), 
        axis.ticks = element_line() )+
  scale_color_manual(values = color_codes) +
  scale_fill_manual( values = color_codes) +
  theme(legend.position = 'n') +
  labs(x = latex2exp::TeX("$\\beta$ in FQ"),
       y = latex2exp::TeX("$\\beta$ in PM")) + 
  geom_text( 
    data = gg1$data[gg1$data$miss>0 & gg1$data$y=='BMI',] %>% 
      {.$th = factor(paste0('missingness > ', .$miss));.}, 
    aes(x = 0.162, y=-0.25, label = paste0("cor = ",round(value,2))), 
    size = 3, color = "black"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 

library(patchwork)
(gg1 + labs(color = "", fill = ""))/gg1s + 
  patchwork::plot_layout(heights = c(1,2.5)) + 
  patchwork::plot_annotation(tag_levels = 'A')

# 13 * 10 
ggsave(filename = 'SFig5_multiplatform_validation_mv.pdf', width = 9, height = 12)


# status
cat("\nFigures saved as PDF in current working directory\n")

