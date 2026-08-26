# Quantile Function of a Reparametrized Distribution

The parent's quantile function, evaluated at the mapped parameters.
Where the parent inverts its distribution function numerically, so does
this; the map costs one evaluation and does not enter the search.

## Usage

``` r
reparam_quantile(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- p:

  A numeric vector of probabilities in \\\[0,1\]\\, or of
  log-probabilities when `log.p` is `TRUE`.

- theta:

  A named list of the new parameters, on the new parameter scale.

- lower.tail:

  Are the probabilities \\P(Y \le q)\\? A single logical, `TRUE` by
  default.

- log.p:

  Is `p` a log-probability? A single logical, `FALSE` by default.

- ...:

  Passed to the parent's method.

## Value

A numeric vector of quantiles, the length of `p` recycled against
`theta`.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md);
[`distrib_cdf.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ReparamContinuousDistrib.md),
which this inverts.

## Examples

``` r
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
distrib_quantile(d, c(0.25, 0.5, 0.75), list(mu = 1, sigma2 = 4))
#> [1] -0.3489795  1.0000000  2.3489795
```
