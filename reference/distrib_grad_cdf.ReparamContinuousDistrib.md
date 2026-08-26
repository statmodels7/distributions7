# Log-CDF Gradient of a Reparametrized Distribution

The chain rule on the parent's cdf derivatives, through
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md).
It is exact whenever the parent's are; where the parent differences its
own cdf, so does this, the chain having nothing closed to carry.

## Arguments

- distrib:

  A `ReparamContinuousDistrib`, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters, on the new parameter scale.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per new parameter, each the length
of `q` recycled against `theta`.

## Details

The map's partials come from the object itself,
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
reading the keyed tables the reparametrization was built with, so a
family created by
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
gets closed cdf derivatives for free as soon as its parent has them.

## See also

[`distrib_hess_cdf.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.ReparamContinuousDistrib.md)
for the second order;
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the identity;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

## Examples

``` r
# A Gaussian reparametrized in the variance. Its parent has closed cdf
# derivatives, so the chain is exact.
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
distrib_grad_cdf(d, c(-1, 0.5, 2), list(mu = 0.3, sigma2 = 1.44))
#> $mu
#> [1] -1.3268963 -0.5790812 -0.1322307
#> 
#> $sigma2
#> [1]  0.59894627 -0.04021397 -0.07805282
#> 
```
