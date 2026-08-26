# Standard Deviations and Correlations of a Multivariate Gaussian

Returns the standard deviations and correlations of the response,
whichever side the parametrization carries, with the closed-form
Jacobian
[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
supplies. A PRECISION parametrization reports two further blocks, which
are what it describes directly: the conditional variances
\\1/\Omega\_{jj} = \operatorname{Var}(Y_j \mid Y\_{-j})\\, and above two
dimensions the partial correlations
\\-\Omega\_{jk}/\sqrt{\Omega\_{jj}\Omega\_{kk}}\\, the correlation of
two coordinates given all the others, which is zero exactly where the
precision has a zero.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- theta:

  A named list of parameters, already aligned by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `value`, `jacobian`, `transform` and `block`, as
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
documents. A covariance parametrization gives \\p(p+1)/2\\ quantities; a
precision one adds \\p\\ conditional variances and, above \\p = 2\\,
\\p(p-1)/2\\ partial correlations.

## What a precision's diagonal means

The quantity with a reading is the conditional VARIANCE, so the diagonal
quantities
[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
produces from \\\Omega\\ are square roots of the wrong thing; they are
dropped and \\1/\Omega\_{jj}\\ is reported instead. Its ratio to the
marginal variance is \\1 - R_j^2\\ for the regression of that coordinate
on all the others.

At \\p = 2\\ there is nothing to condition on, so the partial
correlation IS the correlation and is not printed twice.

## The parametrization's own quantities

Whatever the matrix parametrization declares through
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
is appended as a further block. An AR(1) covariance is about a scale and
a correlation; a log-Cholesky one declares nothing and the summary stops
at the standard deviations.

## Notation

\\\Sigma\\ is the covariance, \\\Omega = \Sigma^{-1}\\ the precision,
\\p\\ the dimension, \\R_j^2\\ the coefficient of determination of the
regression of coordinate \\j\\ on the others, and \\Y\_{-j}\\ the
response with that coordinate removed.

## See also

[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
for the closed-form Jacobian,
[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md)
for the appended block,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the printed result, and
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
mv_derived(d, theta)$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 

# The precision side reports the same law's standard deviations and
# correlation, and adds the conditional variances.
o <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))
th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0, omega_log_L2 = 0,
             omega_L2.1 = 0.5)
od <- mv_derived(o, th_o)
od$value
#>      sd_v1      sd_v2  cor_v1_v2    cvar_v1    cvar_v2 
#>  1.1180340  1.0000000 -0.4472136  1.0000000  0.8000000 
od$block
#>                   sd_v1                   sd_v2               cor_v1_v2 
#>   "Standard deviations"   "Standard deviations"          "Correlations" 
#>                 cvar_v1                 cvar_v2 
#> "Conditional variances" "Conditional variances" 

# A conditional variance is 1 / Omega_jj, and is at most the marginal one.
Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#> conditional    marginal 
#>        1.00        1.25 

# At three dimensions the partial correlations appear as a block of their
# own, the partial and the marginal no longer coinciding.
o3 <- mvgaussian_distrib(3, omega = parameters7::log_cholesky(3))
th3 <- as.list(stats::setNames(
  c(0, 0, 0, 0, 0, 0, 0.5, -0.4, 0.3), o3@params))
unique(mv_derived(o3, th3)$block)
#> [1] "Standard deviations"   "Correlations"          "Conditional variances"
#> [4] "Partial correlations" 
```
