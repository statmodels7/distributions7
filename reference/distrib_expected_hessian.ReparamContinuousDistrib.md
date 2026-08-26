# Expected Information of a Reparametrized Distribution

The same chain rule as the observed Hessian, with the term in the
parent's score dropped: expectation is linear and the map deterministic,
and the score has mean zero, so \$\$E\[\ell^{(ab)}\] = \sum\_{i,j}
E\[\ell^{(ij)}\] \frac{\partial\theta_i}{\partial\psi_a}
\frac{\partial\theta_j}{\partial\psi_b}.\$\$ A parent with an exact
expected information therefore gives an exact one here, and the map's
second partials are not read at all.

## Usage

``` r
reparam_expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, a numeric vector. Its length is what the result is
  recycled to; the values are read only where the parent's own expected
  information reads them.

- theta:

  A named list of the new parameters, on the new parameter scale.

- scale:

  Either `"parameter"` (the default) or `"link"`, handled by the generic
  after this method has returned.

- ...:

  Passed to the parent's method, which is where `approx` and `nsim` are
  read for a family that approximates its expected information.

## Value

A named list of numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `y` recycled against `theta`, in that enumeration's
order.

## Notation

\\\ell\\ is the log-density, \\\theta\\ the parent's parameters and
\\\psi\\ the new ones.

## See also

[`distrib_hessian.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ReparamContinuousDistrib.md),
which keeps the score term;
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md);
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

a <- distrib_expected_hessian(d, c(0, 1, 2), th)
b <- distrib_expected_hessian(gaussian2_distrib(), c(0, 1, 2), th)

# Matched by name, the two enumerating in different orders.
max(abs(unlist(a[names(b)]) - unlist(b)))
#> [1] 1.974032e-12
```
