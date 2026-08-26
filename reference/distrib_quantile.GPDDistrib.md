# Generalized Pareto Quantile Function

Computes the quantiles in closed form, \$\$Q(p; \sigma, \xi) =
\dfrac{\sigma}{\xi}\left\\(1-p)^{-\xi} - 1\right\\,\$\$ with
\\-\sigma\log(1-p)\\ at \\\xi = 0\\. The generalized Pareto is one of
the few families here whose quantile function is elementary, so nothing
is inverted numerically and
[`distrib_rng.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GPDDistrib.md)
can use the inverse transform.

Near zero the shape branch is taken on \\\|\xi\| \< 10^{-8}\\ and the
value comes from `log1p(-p)`, which keeps its digits for a probability
close to zero where `log(1 - p)` would not.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  when `log.p = TRUE`. At `p = 1` the value is the endpoint, which is
  `Inf` for a non-negative shape.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `p`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is the survival probability.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm and is
  exponentiated first. Defaults to `FALSE`.

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

# A negative shape has a finite upper quantile: the endpoint itself.
c(q999 = distrib_quantile(d, 0.999, list(sigma = 2, xi = -0.4)),
  endpoint = distributions7:::gpd_endpoint(2, -0.4))
#>     q999 endpoint 
#> 4.684521 5.000000 
```
