# Fourth-Order Derivatives of a Reparametrized Distribution

Computes the partition sum at order four, observed or expected. A
marginal criterion reads this order in its exact Hessian, so a
reparametrized family used inside a penalized fit reaches it through
this method. The fifth derivative has no method registered on the
reparametrized classes:
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
dispatches to the base method, which differences the lower orders
numerically.

## Usage

``` r
reparam_deriv4(
  distrib,
  y,
  theta,
  expected = FALSE,
  scale = c("parameter", "link"),
  method = NULL,
  ...
)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, a numeric vector.

- theta:

  A named list of the new parameters, on the new parameter scale.

- expected:

  Should the expected derivatives be carried? It is a single logical,
  `FALSE` by default.

- scale:

  Either `"parameter"` (the default) or `"link"`, handled by the generic
  after this method has returned.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the class, so this
  method receives `NULL` and forwards `NULL` to the parent, which
  returns its closed form where it has one. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors keyed as
[`deriv_names(distrib@params, 4)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).
Each has the length of `y` after it is recycled against `theta`.

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the identity;
[`distrib_deriv3.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ReparamContinuousDistrib.md)
for the order below;
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md),
whose accuracy this inherits;
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

# The result has five components for a two-parameter family.
names(distrib_deriv4(d, c(0, 1, 2), th))
#> [1] "mu_mu_mu_mu"                 "mu_mu_mu_sigma2"            
#> [3] "mu_mu_sigma2_sigma2"         "mu_sigma2_sigma2_sigma2"    
#> [5] "sigma2_sigma2_sigma2_sigma2"
```
