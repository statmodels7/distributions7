# Log-CDF Gradient of a Reparametrized Distribution

Applies the chain rule to the parent's cdf derivatives, through
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md).
It is exact whenever the derivatives of the parent and the partials of
the map are; where the parent differences its own cdf, so does this
method, because the chain then carries only numerical derivatives.

## Arguments

- distrib:

  A `ReparamContinuousDistrib`, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters, on the new parameter scale.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per new parameter; each has the
length of `q` recycled against `theta`.

## Details

The partials of the map come from the object itself:
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
reads the keyed tables with which the reparametrization was built. They
are exact when
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
was given `map_derivs`, and single stencils on the map otherwise, so a
family created by
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
carries the closed cdf derivatives of its parent through the map to the
accuracy of those partials.

## See also

[`distrib_hess_cdf.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.ReparamContinuousDistrib.md)
for the second order;
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the identity;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

## Examples

``` r
# This example reparametrizes a Gaussian in the variance. Its parent has
# closed cdf derivatives, and the map's partials here come from a stencil.
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
