# Density of a Reparametrized Distribution

The parent's density, evaluated at the mapped parameters. A change of
parametrization does not change the law, so nothing is recomputed here:
the new parameters go through
[`reparam_theta()`](https://statmodels7.github.io/distributions7/reference/reparam_theta.md)
and the parent answers.

## Usage

``` r
reparam_pdf(distrib, y, theta, log = FALSE, ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  A numeric vector of observations.

- theta:

  A named list of the new parameters, on the new parameter scale.

- log:

  Should the log-density be returned? A single logical, `FALSE` by
  default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, the length of `y` recycled against
`theta`.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md);
[`reparam_theta()`](https://statmodels7.github.io/distributions7/reference/reparam_theta.md)
for the map;
[`distrib_gradient.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ReparamContinuousDistrib.md),
where the parametrization does enter.

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

# The same numbers as the family written out by hand.
all.equal(distrib_pdf(d, c(0, 1, 2), th),
          distrib_pdf(gaussian2_distrib(), c(0, 1, 2), th))
#> [1] TRUE
```
