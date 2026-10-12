# Generalized Pareto Quantile Function

Computes the quantiles in closed form, \$\$Q(p; \sigma, \xi) =
\dfrac{\sigma}{\xi}\left\\(1-p)^{-\xi} - 1\right\\,\$\$ with
\\-\sigma\log(1-p)\\ at \\\xi = 0\\. The quantile function is
elementary, so nothing is inverted numerically and
[`distrib_rng.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GPDDistrib.md)
uses the inverse transform.

The method first forms the logarithm of the survival probability, \\\log
S\\, from `p` without computing \\1 - p\\ or exponentiating: `log1p(-p)`
for a lower-tail probability, \\\log p\\ for an upper-tail one, and the
corresponding forms when `log.p = TRUE`. The quantile is then
\\\sigma\\\mathrm{expm1}(-\xi\log S)/\xi\\, and \\-\sigma\log S\\ at
\\\xi = 0\\, so a probability close to zero or close to one keeps its
digits whichever tail it is given in.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  when `log.p = TRUE`. At `p = 1` the value is the endpoint, which is
  `Inf` for a non-negative shape; a value outside the range gives `NaN`.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `p`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is the survival probability.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm.
  Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of the length of the recycled inputs.

## See also

[`distrib_cdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GPDDistrib.md),
which it inverts,
[`distrib_rng.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GPDDistrib.md),
which draws from it, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
th <- list(sigma = 1.5, xi = 0.3)

# The round trip is exact, both functions being closed form.
y <- c(0.2, 1, 4)
all.equal(distrib_quantile(d, distrib_cdf(d, y, th), th), y)
#> [1] TRUE

# Shape zero is the exponential quantile.
p <- c(0.1, 0.5, 0.9)
all.equal(distrib_quantile(d, p, list(sigma = 1.5, xi = 0)),
          qexp(p, rate = 1 / 1.5))
#> [1] TRUE

# A negative shape bounds the quantiles by the endpoint, which p = 1
# reaches.
c(q999 = distrib_quantile(d, 0.999, list(sigma = 2, xi = -0.4)),
  q1 = distrib_quantile(d, 1, list(sigma = 2, xi = -0.4)),
  endpoint = distributions7:::gpd_endpoint(2, -0.4))
#>     q999       q1 endpoint 
#> 4.684521 5.000000 5.000000 
```
