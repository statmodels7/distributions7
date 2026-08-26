# Moments of a Reparametrized Distribution

[`mean()`](https://rdrr.io/r/base/mean.html),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
all delegate to the parent at the mapped parameters. A reparametrization
does not change the law, so it does not change a moment; what changes is
the coordinates the moment is computed from. Delegating keeps the
parent's closed forms, where falling through to the `distrib` default
would run a quadrature.

## Usage

``` r
reparam_mean(x, theta, ...)

reparam_variance(x, theta, ...)

reparam_skewness(x, theta, ...)

reparam_kurtosis(x, theta, ...)
```

## Arguments

- x:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- theta:

  A named list of the new parameters, on the new parameter scale.

- ...:

  Passed to the parent's method, and from there to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  if the parent has no closed form.

## Value

A numeric vector: means for
[`mean()`](https://rdrr.io/r/base/mean.html), variances for
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
and the standardized third and fourth moments for
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md),
the excess in the last case. `NaN` or `Inf` wherever the parent's own
method gives one.

## See also

[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the numerical route the parent may take;
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
for the default this replaces;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

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
th <- list(mu = 1, sigma2 = 4)

# The parent's closed forms, in the new coordinates.
c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#> [1] 1 4 0 0
```
