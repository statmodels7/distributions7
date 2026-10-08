# Generalized Pareto Log-CDF Derivatives

Return the derivatives in \\(\sigma, \xi)\\ of the distribution function
or of the survival function, or of their logarithms, at orders one to
four, each from its own compiled kernel.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value), each a numeric vector of length 1 or of the length of `q`.

- lower.tail:

  Logical; if `TRUE` (the default) the derivatives of \\F\\, otherwise
  those of \\S = 1 - F\\.

- log:

  Logical; if `TRUE` (the default) the derivatives of the logarithm.

- ...:

  Unused.

## Value

A named list of numeric vectors keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md):
two components for the gradient, three for the Hessian, four at order
three and five at order four.

## Details

The survival function is \\S = e^{-W}\\ with \\W = z\\\phi(\xi z)\\ as
[`distrib_deriv3.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
writes it, so \\\log S = -W\\ and its derivatives are those of \\W\\,
free of any division by \\\xi\\. The derivatives of \\\log F = \log(1 -
e^{-W})\\ follow by Faa di Bruno with \\G(w) = \log(1 - e^{-w})\\, whose
derivatives are polynomials in \\q = 1/(e^{W} - 1)\\: \\G' = q\\ and
\\dq/dw = -q(1+q)\\. Those of \\S\\ itself are \\S\\ times the
corresponding polynomial in the derivatives of \\-W\\, and \\\partial F
= -\partial S\\.

At \\\xi \ge 0\\ the support is \\(0, \infty)\\; at \\\xi \< 0\\ it is
bounded above at \\\sigma/\|\xi\|\\. Past the upper end, and below zero,
\\F\\ is constant and every derivative is zero, except the derivatives
of \\\log F\\ at \\q \le 0\\, where \\\log F = -\infty\\ and they are
`NaN`.

## See also

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
