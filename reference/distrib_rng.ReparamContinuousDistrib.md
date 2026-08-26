# Random Generation from a Reparametrized Distribution

The parent's generator, at the mapped parameters. The draws come from
the parent's own method and consume its random numbers, so a seed set
before the call gives the same sample the parent would give at the
mapped values.

## Usage

``` r
reparam_rng(distrib, n, theta, ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- n:

  The number of draws, a single non-negative whole number.

- theta:

  A named list of the new parameters, on the new parameter scale. A
  component of length `n` gives one draw per setting.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md);
[`distrib_quantile.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ReparamContinuousDistrib.md),
which an inverse-transform parent uses.

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

# The same draws the parent gives at the mapped parameters.
set.seed(1); a <- distrib_rng(d, 5, list(mu = 1, sigma2 = 4))
set.seed(1); b <- distrib_rng(gaussian1_distrib(), 5,
                              list(mu = 1, sigma = 2))
all.equal(a, b)
#> [1] TRUE
```
