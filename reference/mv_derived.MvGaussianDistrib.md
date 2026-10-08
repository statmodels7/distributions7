# Interpretable Quantities of a Multivariate Gaussian

Returns the quantities the chosen parametrization describes, with their
Jacobian. Where the matrix parametrization declares its own block
through
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
(a compound symmetry, an AR(1), an AR(p)), that block alone is returned,
with the common variance `scale` reported as its square root `sd` on the
covariance side. Otherwise a COVARIANCE parametrization reports the
standard deviations and the correlations, and a PRECISION
parametrization reports the conditional variances \\1/\Omega\_{jj} =
\operatorname{Var}(Y_j \mid Y\_{-j})\\ and the partial correlations
\\-\Omega\_{jk}/\sqrt{\Omega\_{jj}\Omega\_{kk}}\\, the correlation of
two coordinates given all the others.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
  or
  [`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- theta:

  A named list of parameters, already aligned by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `value`, `jacobian`, `transform` and `block`, as
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
documents: the structure's own block where it declares one; otherwise
\\p\\ standard deviations and \\p(p-1)/2\\ correlations for a
covariance, or \\p\\ conditional variances and \\p(p-1)/2\\ partial
correlations for a precision.

## One parametrization, one set of quantities

A user who writes the model on the precision asks to read the precision,
so the marginal standard deviations and correlations are not reported
beside its own quantities; they remain available from
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md).
In the same way a structured matrix is fixed by its few parameters, and
the standard deviation of every coordinate and the correlation of every
pair would repeat them.

## What a precision's diagonal means

The quantity with a reading is the conditional VARIANCE, so
\\1/\Omega\_{jj}\\ is reported rather than a square root of
\\\Omega\_{jj}\\. Its ratio to the marginal variance is \\1 - R_j^2\\
for the regression of that coordinate on all the others. At \\p = 2\\
the partial correlation equals the correlation.

## Notation

\\\Sigma\\ is the covariance, \\\Omega = \Sigma^{-1}\\ the precision,
\\p\\ the dimension, \\R_j^2\\ the coefficient of determination of the
regression of coordinate \\j\\ on the others, and \\Y\_{-j}\\ the
response with that coordinate removed.

## See also

[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
for the closed-form Jacobian,
[`mv_own_block()`](https://statmodels7.github.io/distributions7/reference/mv_own_block.md)
for the structure's own block,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the printed result, and
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
mv_derived(d, theta)$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 

# The precision side reports the conditional variances and the partial
# correlation, and nothing of the covariance.
o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0, omega_log_L2 = 0,
             omega_L2.1 = 0.5)
od <- mv_derived(o, th_o)
od$value
#>    cvar_v1    cvar_v2 pcor_v1_v2 
#>  1.0000000  0.8000000 -0.4472136 

# A conditional variance is 1 / Omega_jj, and is at most the marginal one.
Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#> conditional    marginal 
#>        1.00        1.25 

# A compound symmetry reports its standard deviation and its correlation.
cs <- mvgaussian1_distrib(3, parameters7::compound_symmetry(3))
mv_derived(cs, list(mu1 = 0, mu2 = 0, mu3 = 0, sigma_log_scale = log(4),
                    sigma_logit_rho = 0))$value
#>   sd  rho 
#> 2.00 0.25 
```
