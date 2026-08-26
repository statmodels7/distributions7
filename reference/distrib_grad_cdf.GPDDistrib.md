# Generalized Pareto Log-CDF Derivatives

Closed form at every order from one to four, from the survival function
\\S = (1 + \xi q/\sigma)^{-1/\xi}\\. Its logarithm is written \\L =
-(q/\sigma)\\\Lambda(\xi q/\sigma)\\ with \\\Lambda(u) = \log(1+u)/u\\,
which carries no division by the shape, so the exponential limit \\\xi
\to 0\\ is an ordinary point of the formula and not a branch.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- q:

  A numeric vector of quantiles. Values outside the support give
  derivatives of exactly zero.

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of the order the generic asked for,
keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md):
two components for the gradient, three for the Hessian, four at order 3
and five at order 4.

## The support moves with the shape

At \\\xi \ge 0\\ the support is \\(0, \infty)\\; at \\\xi \< 0\\ it is
bounded above at \\\sigma/\|\xi\|\\, and past that endpoint every
derivative is exactly zero. The mask is computed in
[`gpd_surv_pieces()`](https://statmodels7.github.io/distributions7/reference/gpd_surv_pieces.md),
the family's fixed bounds being unable to record a support that moves
with a parameter.

## What it is worth, and the limit as a check

Against a product stencil on the same cdf at \\\sigma = 1\\, \\\xi =
0.3\\: \\5.2\times10^{-11}\\ at order 1, \\2.7\times10^{-7}\\ at order
2, \\1.7\times10^{-5}\\ at order 3 and \\8.4\times10^{-4}\\ at order 4.
At \\\xi = 0\\ the scale component equals the exponential family's to
the last bit, and the fourth derivative reads the same value at \\\xi =
10^{-7}\\ and at \\10^{-9}\\, which is what a removable singularity
handled properly looks like.

## Notation

\\\sigma \> 0\\ is the scale, \\\xi\\ the shape of either sign, \\u =
\xi q/\sigma\\, \\\Lambda(u) = \log(1+u)/u\\, \\F\\ the distribution
function and \\S = 1 - F\\ the survival function.

## See also

[`gpd_surv_pieces()`](https://statmodels7.github.io/distributions7/reference/gpd_surv_pieces.md)
and
[`gpd_lambda_derivs()`](https://statmodels7.github.io/distributions7/reference/gpd_lambda_derivs.md)
for the construction;
[`distrib_grad_cdf.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.ExponentialDistrib.md),
the \\\xi = 0\\ case;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## Examples

``` r
d <- gpd_distrib()
q <- c(0.5, 2, 5)

# At shape zero the scale component is the exponential family's.
rbind(gpd = distrib_grad_cdf(d, q, list(sigma = 3, xi = 0),
                             log = FALSE)$sigma,
      exponential = distrib_grad_cdf(exponential_distrib(), q,
                                     list(mu = 3), log = FALSE)$mu)
#>                    [,1]       [,2]       [,3]
#> gpd         -0.04702676 -0.1140927 -0.1049309
#> exponential -0.04702676 -0.1140927 -0.1049309

# A negative shape bounds the support at sigma / |xi| = 2.
distrib_grad_cdf(d, c(1, 2, 3), list(sigma = 1, xi = -0.5),
                 log = FALSE)$sigma
#> [1] -0.5  0.0  0.0
```
