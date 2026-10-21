# Student t Log-CDF Hessian

Computes the Hessian in closed form in the location-scale block, and in
the mixed components with a shape \\k\\, which are \\-f s_k\\ and \\-z f
s_k\\ at \\q\\; the components that involve only shape parameters are
integrals of the density's own derivatives, taken by the compiled rule.
The method is
[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
itself.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive); each is a numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

## Value

A named list of six numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md);
each has the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
degrees of freedom, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
for the shared body;
[`distrib_grad_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.StudentT1Distrib.md)
for the first order;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()
th <- list(mu = 0.3, sigma = 1.2, nu = 6)
q <- c(-1, 0.5, 2)

# The result has six components, keyed as hess_names().
names(distrib_hess_cdf(d, q, th))
#> [1] "mu_mu"       "sigma_sigma" "nu_nu"       "mu_sigma"    "mu_nu"      
#> [6] "sigma_nu"   
```
