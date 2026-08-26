# Pseudo-Huber Quantile Function

Computes \\Q(p)\\ by root-finding on the numerical distribution
function, the family having no elementary quantile. Symmetry about
\\\mu\\ does most of the work: \\Q(1/2) = \mu\\ exactly and without a
search, and a probability above one half is reflected through \\Q(p) =
2\mu - Q(1-p)\\, so the search always runs in the lower half.

The bracket starts at ten standard deviations below the location and
doubles its width until it contains the root, at most 100 times.
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) then closes
it to `sqrt(.Machine$double.eps)`. Each evaluation of the objective is a
quadrature, so this is by far the dearest method of the family.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. `NaN` is returned for `NA` and for a
  value outside the range; 0 gives `-Inf` and 1 gives `Inf`.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `p`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is read as a logarithm. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles on \\\[-\infty, \infty\]\\, of length
`max(length(p), length(mu), length(sigma), length(nu))`.

## See also

[`distrib_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuberDistrib.md)
for the function inverted here,
[`distrib_rng.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuberDistrib.md),
which draws by inverting it at uniform variates, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The median is the location, returned without a search.
distrib_quantile(d, 0.5, th)
#> [1] 0.4

# The quartiles, and the round trip back through the distribution function.
q <- distrib_quantile(d, c(0.25, 0.75), th)
q
#> [1] -0.8235709  1.6235709
distrib_cdf(d, q, th)
#> [1] 0.25 0.75

# Symmetry: the two quartiles are equidistant from the location.
q - 0.4
#> [1] -1.223571  1.223571
```
