# Chi-Squared Fifth-Order Derivatives

Computes the fifth derivative of the chi-squared log-density with
respect to its degrees of freedom, in closed form from a compiled
kernel, \$\$\ell^{(\mu\mu\mu\mu\mu)} =
-\tfrac{1}{32}\\\psi^{(4)}(\mu/2),\$\$ which does not depend on the
response. With `scale = "link"` the components are composed with the
link derivatives by the routine that every family with a fifth-order
kernel uses.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`.

- scale:

  One of `"parameter"` (the default) or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list with one numeric vector, `mu_mu_mu_mu_mu`, of length
`max(length(y), length(mu))`.

## See also

[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
all.equal(distrib_deriv5(d, c(0.5, 2, 6), list(mu = 3))$mu_mu_mu_mu_mu,
          rep(-psigamma(1.5, 4) / 32, 3))
#> [1] TRUE
```
