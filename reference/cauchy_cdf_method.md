# Cauchy Distribution Function Derivatives

Returns the method body that
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
and
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
use for the Cauchy, at one order. With \\Y = q - \mu + i\sigma\\, \\F =
1 - \mathrm{Im}(\log Y)/\pi\\, so every derivative of order \\n\\ with
\\j\\ derivatives in \\\mu\\ and \\k\\ in \\\sigma\\ is one term,
\$\$-\tfrac{1}{\pi}(-1)^{j+n-1}(n-1)!\\\mathrm{Im}(i^k Y^{-n}),\$\$
computed by a compiled kernel without squaring the residual. On the log
scale the method divides by the probability of the requested tail, which
[`stats::pcauchy()`](https://rdrr.io/r/stats/Cauchy.html) returns for
either tail.

## Usage

``` r
cauchy_cdf_method(order)
```

## Arguments

- order:

  The derivative order, 1 to 4.

## Value

A function with the signature of the four generics, returning a named
list keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
or as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2.

## Examples

``` r
d <- cauchy_distrib()
th <- list(mu = 0, sigma = 1)
# Far in the tail the log-scale derivatives in mu tend to those of
# log(1/(pi q)).
distrib_grad_cdf(d, 1e8, th, lower.tail = FALSE)$mu
#> [1] 1e-08
```
