# Log-Survival Derivatives of the Weibull, Mean and Shape

Computes the derivatives of order `o` of \\\log S(q) = -e^{h}\\, \\h =
\sigma\\\log q - \log m + \log\Gamma(1 + 1/\sigma)\\\\, in \\(m,
\sigma)\\, as \\-e^{h}\\ times the complete Bell polynomial in the
partial derivatives of \\h\\
([`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)).
The function \\h\\ is \\\sigma\\ times a sum of a function of \\m\\ and
a function of \\\sigma\\, so a mixed partial derivative of \\h\\ is
nonzero only when it names \\\sigma\\ at most once; the pure derivatives
in \\\sigma\\ are \\\sigma g_k + k g\_{k-1}\\, where \\g_k\\ is the
derivative of order \\k\\ of \\\log\Gamma(1 + 1/\sigma)\\ by the chain
rule over the polygamma functions and \\g_0\\ is the whole sum. At \\q
\le 0\\ every derivative is zero.

## Usage

``` r
weibull3_log_surv_derivs(distrib, q, theta, o)
```

## Arguments

- distrib:

  A `Weibull3Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mean` and `sigma`.

- o:

  The derivative order, 1 to 4.

## Value

A named list, one component per multi-index of order `o`, in the order
of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order two.

## Examples

``` r
d <- weibull3_distrib()
distributions7:::weibull3_log_surv_derivs(d, c(18, 480),
                                          list(mean = 60, sigma = 6), 2L)
#> $mean_mean
#> [1] -5.422180e-06 -1.949783e+03
#> 
#> $sigma_sigma
#> [1] -6.987221e-04 -7.100705e+05
#> 
#> $mean_sigma
#> [1] -4.912391e-05  3.720917e+04
#> 
```
