# Is This a Reparametrized Distribution?

`TRUE` for either of the two wrapper classes,
[`ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md)
and
[`ReparamDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md),
and `FALSE` for anything else. The test exists because the two classes
have different parents, so no single `S7_inherits()` call answers it.

## Usage

``` r
is_reparam(distrib)
```

## Arguments

- distrib:

  An object inheriting from `distrib`, or anything else.

## Value

A single logical.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
which builds them;
[`ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md)
for the classes.

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
c(distributions7:::is_reparam(d),
  distributions7:::is_reparam(gaussian2_distrib()))
#> [1]  TRUE FALSE
```
