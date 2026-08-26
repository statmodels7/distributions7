# Log-CDF Hessian of a Reparametrized Distribution

The second-order chain rule on the parent's cdf derivatives, through
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md):
\\\partial^2 F/\partial\psi_a\partial\psi_b = \sum\_{k,l} F\_{kl} h^k_a
h^l_b + \sum_k F_k h^k\_{ab}\\. The second sum is the one a first-order
chain does not have, and it is why the map's second partials are needed
as well as its first.

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

A named list of numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`. The gradient is not
returned alongside.

## Details

The gate asks the parent for exactness at **both** orders. A parent with
a closed gradient and a differenced Hessian sends this whole method to
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md),
because carrying an exact first order onto an approximate second would
not make the second exact.

## See also

[`distrib_grad_cdf.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.ReparamContinuousDistrib.md)
for the first order;
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the identity;
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
q <- c(-1, 0.5, 2)
th <- list(mu = 0.3, sigma2 = 1.44)

# Against a central difference of the reparametrized cdf.
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) - unlist(fd)))
#> [1] 6.091732e-09
```
