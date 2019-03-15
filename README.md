
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/cmpreg_hex.png" align="right" height=160/ display="block">

# `cmpreg`: Reparametrized COM-Poisson Regression Models

[![Travis build
status](https://travis-ci.org/jreduardo/cmpreg.svg?branch=master)](https://travis-ci.org/jreduardo/cmpreg)

> [Eduardo E. R. Junior](http://leg.ufpr.br/~eduardojr) -
> <jreduardo@usp.br>, IME-USP

The `cmpreg` package contains functions to fit Conway-Maxwell-Poisson
(COM-Poisson) models with varying dispersion (model mean and dispersion
parameters as functions of covariates.) in the mean-type parametrization
proposed by Ribeiro Jr, et al. 2018 \<arxiv.org/abs/1801.09795\>. The
functions to computate the normalizing constant are written in C++, so
the computations is reasonably fast.

Joint work with [Walmes M. Zeviani](www.leg.ufpr.br/~walmes/) and
[Clarice G.B.
Demétrio](http://ce.esalq.usp.br/equipe/clarice-garcia-borges-demetrio).

## Installation

You can install the development version of `cmpreg` from
[GitHub](https://github.com/jreduardo/cmpreg) with:

``` r

# install.packages("devtools")
devtools::install_github("jreduardo/cmpreg")
```

## Usage and example

Basically, this package implements methods similar to those related to
glm objects. The main function is `cmp(...)`.

``` r

library(cmpreg)

# Fit model ------------------------------------------------------------
model <- cmp(formula = ninsect ~ extract,
             dformula = ~extract,
             data = sitophilus)

# Methods --------------------------------------------------------------
print(model)
#> 
#> COM-Poisson regression models
#> Call:  cmp(formula = ninsect ~ extract,
#>            dformula = ~extract,
#>            data = sitophilus)
#> 
#> Mean coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>      3.449860      -0.006594      -0.052379      -3.310863  
#> 
#> Dispersion coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>       -0.6652        -0.3831        -0.3724        -0.1186  
#> 
#> Residual degrees of freedom: 32
#> Minus twice the log-likelihood: 242.8278
summary(model)
#> 
#> Individual Wald-tests for COM-Poisson regression models
#> Call:  cmp(formula = ninsect ~ extract,
#>            dformula = ~extract,
#>            data = sitophilus)
#> 
#> Mean coefficients:
#>                Estimate Std. Error Z value Pr(>|z|)    
#> (Intercept)    3.449860   0.077995  44.232  < 2e-16 ***
#> extractLeaf   -0.006594   0.122209  -0.054    0.957    
#> extractBranch -0.052379   0.123463  -0.424    0.671    
#> extractSeed   -3.310863   0.543822  -6.088 1.14e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Dispersion coefficients:
#>               Estimate Std. Error Z value Pr(>|z|)
#> (Intercept)    -0.6652     0.4573  -1.455    0.146
#> extractLeaf    -0.3831     0.6509  -0.589    0.556
#> extractBranch  -0.3724     0.6514  -0.572    0.567
#> extractSeed    -0.1186     1.5502  -0.077    0.939
#> 
#> Residual degrees of freedom: 32
#> Minus twice the log-likelihood: 242.8278
equitest(model)
#> 
#> Likelihood ratio test for equidispersion 
#> 
#>         Resid.df   Loglik LRT_stat LRT_df Pr(>LRT_stat)   
#> Model 1       32 -121.414   18.338      4       0.00106 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Predict new data -----------------------------------------------------
newdf <- sitophilus[c(1, 11, 21, 31), -2, drop = FALSE]
predict(model,
        newdata = newdf,
        what = "all",
        type = "response",
        se.fit = TRUE,
        augment_data = TRUE)
#>   extract       what        fit       ste
#> 1    Leaf       mean 31.2889797 2.9437657
#> 2  Branch       mean 29.8887063 2.8605423
#> 3    Seed       mean  1.1491206 0.3120492
#> 4 Control       mean 31.4959708 2.4565258
#> 5    Leaf dispersion  0.3505199 0.1623447
#> 6  Branch dispersion  0.3542934 0.1643381
#> 7    Seed dispersion  0.4566733 0.3412895
#> 8 Control dispersion  0.5141727 0.2351421
```

Currently, the methods implemented for `"cmpreg"` objects are

``` r

methods(class = "cmpreg")
#>  [1] anova        coef         equitest     fitted       logLik      
#>  [6] model.matrix predict      print        summary      vcov        
#> see '?methods' for accessing help and source code
```

## Related projects

There are other R packages to deal with COM-Poisson models that have
somehow contributed to the writing of `cmpreg`.

  - distribution under original parametrization.
  - [`CompGLM`](https://github.com/jeffpollock9/CompGLM): Fit
    COM-Poisson models under original parametrization (includes
    dispersion modeling).
  - parametrization (includes zero-inflation and dispersion modeling).
  - [`glmmTMB`](https://github.com/glmmTMB/glmmTMB): Fit (among other)
    COM-Poisson models under a different mean-parametrization (includes
    zero-inflation, dispersion modeling and random effects).

## License

The `gammacount` package is licensed under the [GNU General Public
License, version 3](https://www.gnu.org/licenses/gpl-3.0.ht), descrita
no arquivo `LICENSE.md`, © 2019 E. E., Ribeiro Jr.

<!------------------------------------------- -->

<!-- Links -->
