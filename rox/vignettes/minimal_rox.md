# Minimal rox example

This tutorial outlines the functionality of `rox`.

``` r
# initialize
set.seed(42)
library(survival)
library(rox)

# background population size
N = 1e4

# let y be an outcome of interest
y = rnorm(N)

# let x be a metabolite measurement
x = y + rnorm(N)

# ground truth: concordance between x and y
d0 = survival::concordance(x~y) # we are using the concordance function from the survival package
d0$concordance
```

    ## [1] 0.7489285

``` r
# limit-of-detection (LOD) for missing values
th = 0

# n: size of sample from the larger population for analysis
n = 400

# draw sample
i = sample(N, n)
xi = x[i]
yi = y[i]

# set values below th to missing
xh = xi
xh[xh<th] = NA

# show data
data.table::data.table(xh = xh, yi)
```

    ##              xh          yi
    ##   1: 1.16215361 -0.98713243
    ##   2: 0.59026876  0.11335259
    ##   3:         NA -0.85343702
    ##   4: 0.34359125 -0.33236346
    ##   5: 3.69386715  1.87516638
    ##  ---                       
    ## 396:         NA -2.32709866
    ## 397: 0.39789320  0.02847577
    ## 398: 0.04149851  2.03488046
    ## 399: 0.28997557 -0.97186164
    ## 400: 0.13566929  0.76186345

``` r
# rox: estimates the concordance with missing values
rox(xh, yi)
```

    ## 
    ##      d    se      z     p
    ##   0.7431 0.0188 12.9018 <2e-16
    ## ---
    ## w = 0.5

``` r
# sample estimate if there were no missing values (should be close to ground truth d0)
(d1 = concordance(xi~yi))
```

    ## Call:
    ## concordance.formula(object = xi ~ yi)
    ## 
    ## n= 400 
    ## Concordance= 0.7433 se= 0.01289
    ## concordant discordant     tied.x     tied.y    tied.xy 
    ##      59317      20483          0          0          0

``` r
# compare
rbind(
  `concordance with no missing values:` = unname(d1$concordance),
  `rox estimate with missing values:` = unname(rox(xh, yi)$stats["d"])
)
```

    ##                                          [,1]
    ## concordance with no missing values: 0.7433208
    ## rox estimate with missing values:   0.7431328



## Multivariable example

``` r

# population size
N = 1e4 # background
n = 500 # analyzed sample

# let y be an outcome of interest
y = rnorm(N)

# let z be a covariate of interest
z = rnorm(N)

# let x be a metabolite measurement
x = y + 2*z + rnorm(N)

# fit for background population = ground truth
fit0 = rox_mv(x, data.frame(y,z))
fit0
```

    ## 
    ## (coefficients for multivariable model) 
    ## model: concreg 
    ##       B      se      z p
    ## y 1.223 0.01615  75.72 0
    ## z 2.513 0.02432 103.34 0
    ## ---------------------------------------
    ## 
    ##      d    se      z     p
    ##   0.8662 0.0026 142.9032 <2e-16
    ## ---
    ## no missing value.

``` r

# draw a sample from the population
i = sample(N, n)
xi = x[i]
yi = y[i]
zi = z[i]

# set values below th to missing
xh = xi
xh[xh<th] = NA

# show data
data.table::data.table(xh = xh, yi, zi)
```

    ##             xh         yi         zi
    ##   1:        NA  0.2497746  0.3372729
    ##   2:        NA -1.8390140 -0.2302251
    ##   3: 0.8680995 -0.8864904  1.8970968
    ##   4: 3.4546689  1.6819432  0.6521481
    ##   5:        NA -1.1845607 -0.4987884
    ##  ---                                
    ## 496:        NA -0.7834889 -1.1972024
    ## 497: 0.2429925  1.3405111 -0.1658415
    ## 498: 2.1589406 -0.1519269  0.9440151
    ## 499:        NA -0.9620413  0.8673903
    ## 500: 1.0091711  1.7540953 -0.4469757

``` r
# rox: estimates the concordance with missing values
rox_mv(xh, data.frame(y=yi, z=zi))
```

    ## 
    ## (coefficients for multivariable model) 
    ## model: concreg 
    ##       B      se     z         p
    ## y 1.185 0.09219 12.85 8.637e-38
    ## z 2.437 0.14569 16.73 8.037e-63
    ## ---------------------------------------
    ## 
    ##      d    se      z     p
    ##   0.8589 0.0147 24.366 <2e-16
    ## ---
    ## w = 0.522

``` r
# concordance regression estimate with no missing values
(fit1 = rox_mv(xi, data.frame(y=yi, z=zi)))
```

    ## 
    ## (coefficients for multivariable model) 
    ## model: concreg 
    ##       B      se     z          p
    ## y 1.216 0.06932 17.54  7.412e-69
    ## z 2.405 0.10823 22.22 2.267e-109
    ## ---------------------------------------
    ## 
    ##      d    se      z     p
    ##   0.8613 0.0115 31.3513 <2e-16
    ## ---
    ## no missing value.

``` r
# compare betas
rbind(`concordance regression with no missing values:` =
        fit1$stats_model[,"B"],
      `rox coefficient estimates with missing values:` =
        rox_mv(xh, data.frame(y=yi, z=zi))$stats_model[,"B"])
```

    ##                                                       y        z
    ## concordance regression with no missing values: 1.215617 2.404662
    ## rox coefficient estimates with missing values: 1.184582 2.437202
