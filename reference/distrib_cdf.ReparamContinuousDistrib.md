# Distribution Function of a Reparametrized Distribution

The parent's distribution function, evaluated at the mapped parameters.
The `lower.tail` and `log.p` arguments are passed through and applied by
the parent, so they are not applied twice.

## Usage

``` r
reparam_cdf(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters, on the new parameter scale.

- lower.tail:

  Should probabilities be \\P(Y \le q)\\? A single logical, `TRUE` by
  default.

- log.p:

  Should log-probabilities be returned? A single logical, `FALSE` by
  default.

- ...:

  Passed to the parent's method.

## Value

A numeric vector of probabilities, the length of `q` recycled against
`theta`.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md);
[`distrib_grad_cdf.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.ReparamContinuousDistrib.md)
for its derivatives.

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
distrib_cdf(d, c(0, 1, 2), list(mu = 1, sigma2 = 4))
#> [1] 0.3085375 0.5000000 0.6914625
```
