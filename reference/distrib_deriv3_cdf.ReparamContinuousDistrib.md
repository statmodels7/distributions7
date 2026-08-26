# Third and Fourth Log-CDF Derivatives of a Reparametrized Distribution

The chain rule on the parent's cdf derivatives at orders 3 and 4,
through
[`mapped_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv_k.md).
It is exact whenever the parent's are exact at every order up to the one
wanted, and falls to one product stencil on the reparametrized cdf
otherwise.

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

For
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md),
a named list of third-derivative components keyed as
[`deriv_names(distrib@params, 3)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md);
for
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md),
the fourth-order components. Each vector is the length of `q` recycled
against `theta`.

## Details

The map's partials come from the object itself,
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
reading the keyed tables the reparametrization was built with. A family
created by
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
over a parent with closed cdf derivatives therefore reaches the fourth
order with no arithmetic of its own.

## See also

[`chain_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv_k.md)
for the identity;
[`distrib_grad_cdf.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.ReparamContinuousDistrib.md)
for the first order;
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
th <- list(mu = 0.3, sigma2 = 1.44)

# It agrees with the family written out by hand, to the closed route's own
# accuracy rather than to the stencil's.
a <- distrib_deriv3_cdf(d, 1, th)
b <- distrib_deriv3_cdf(gaussian2_distrib(), 1, th)
max(abs(unlist(a) - unlist(b)))
#> [1] 3.766547e-08
```
