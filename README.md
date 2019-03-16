
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
#>      3.449861      -0.006596      -0.052377      -3.311192  
#> 
#> Dispersion coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>       -0.6652        -0.3832        -0.3724        -0.1177  
#> 
#> Residual degrees of freedom: 32
#> Minus twice the log-likelihood: 242.8279

summary(model)
#> 
#> Individual Wald-tests for COM-Poisson regression models
#> Call:  cmp(formula = ninsect ~ extract,
#>            dformula = ~extract,
#>            data = sitophilus)
#> 
#> Mean coefficients:
#>                Estimate Std. Error Z value Pr(>|z|)    
#> (Intercept)    3.449861   0.077995  44.232  < 2e-16 ***
#> extractLeaf   -0.006596   0.122210  -0.054    0.957    
#> extractBranch -0.052377   0.123462  -0.424    0.671    
#> extractSeed   -3.311192   0.541399  -6.116  9.6e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Dispersion coefficients:
#>               Estimate Std. Error Z value Pr(>|z|)
#> (Intercept)    -0.6652     0.4573  -1.455    0.146
#> extractLeaf    -0.3832     0.6509  -0.589    0.556
#> extractBranch  -0.3724     0.6514  -0.572    0.568
#> extractSeed    -0.1177     1.5464  -0.076    0.939
#> 
#> Residual degrees of freedom: 32
#> Minus twice the log-likelihood: 242.8279

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
#> 1    Leaf       mean 31.2889190 2.9438074
#> 2  Branch       mean 29.8887880 2.8605146
#> 3    Seed       mean  1.1487432 0.3120589
#> 4 Control       mean 31.4959985 2.4565378
#> 5    Leaf dispersion  0.3505090 0.1623430
#> 6  Branch dispersion  0.3542998 0.1643400
#> 7    Seed dispersion  0.4570660 0.3423450
#> 8 Control dispersion  0.5141684 0.2351420
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

  - [`compoisson`](https://github.com/cran/compoisson): Routines for
    density and moments of the COM-Poisson distribution under original
    parametrization.
  - [`CompGLM`](https://github.com/jeffpollock9/CompGLM): Fit
    COM-Poisson models under original parametrization (includes
    dispersion modeling).
  - [`COMPoissonReg`](https://github.com/lotze/COMPoissonReg): Fit
    COM-Poisson models under original parametrization (includes
    zero-inflation and dispersion modeling).
  - [`glmmTMB`](https://github.com/glmmTMB/glmmTMB): Fit (among other)
    COM-Poisson models under a different mean-parametrization (includes
    zero-inflation, dispersion modeling and random effects).

## License

The `gammacount` package is licensed under the [GNU General Public
License, version 3](https://www.gnu.org/licenses/gpl-3.0.ht), see file
`LICENSE.md`, © 2019 E. E., Ribeiro Jr.

<!------------------------------------------- -->

<!-- Links -->
