# Observed Hessian of a Reparametrized Distribution

The second-order chain rule, which keeps the term in the parent's score
as well as the one in its Hessian: \$\$\ell^{(ab)} = \sum\_{i,j}
\ell^{(ij)} \frac{\partial\theta_i}{\partial\psi_a}
\frac{\partial\theta_j}{\partial\psi_b} + \sum_i \ell^{(i)}
\frac{\partial^2\theta_i}{\partial\psi_a \partial\psi_b}.\$\$ The second
sum is what a first-order chain does not have, and it is why the map's
second partials are needed as well as its first.

## Usage

``` r
reparam_hessian(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, a numeric vector.

- theta:

  A named list of the new parameters, on the new parameter scale.

- scale:

  Either `"parameter"` (the default) or `"link"`, handled by the generic
  after this method has returned.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `y` recycled against `theta`. Note the **order**:
this route enumerates as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
does, diagonal first, while a family written out by hand may enumerate
lexicographically. Compare two such results by name and not by position.

## Notation

\\\ell\\ is the log-density, \\\theta\\ the parent's parameters and
\\\psi\\ the new ones.

## See also

[`distrib_gradient.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ReparamContinuousDistrib.md)
for the first order;
[`distrib_expected_hessian.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ReparamContinuousDistrib.md),
where the second sum drops;
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md).

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

a <- distrib_hessian(d, c(0, 1, 2), th)
b <- distrib_hessian(gaussian2_distrib(), c(0, 1, 2), th)

# The two enumerate their components in different orders, so match by name.
rbind(reparametrized = names(a), written_out = names(b))
#>                [,1]    [,2]            [,3]           
#> reparametrized "mu_mu" "sigma2_sigma2" "mu_sigma2"    
#> written_out    "mu_mu" "mu_sigma2"     "sigma2_sigma2"
max(abs(unlist(a[names(b)]) - unlist(b)))
#> [1] 1.974032e-12
```
