# Student t Log-CDF Hessian

Closed form in the location-scale block, `mu_mu`, `sigma_sigma` and
`mu_sigma`; the three components touching the degrees of freedom are
differenced, that derivative having no elementary form. The method is
[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
itself, shared with the pseudo-Huber and the skew t.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

## Value

A named list of six numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
degrees of freedom, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
for the shared body;
[`distrib_grad_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.StudentT1Distrib.md)
for the first order;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()
th <- list(mu = 0.3, sigma = 1.2, nu = 6)
q <- c(-1, 0.5, 2)

# Six components: three closed, three differenced.
names(distrib_hess_cdf(d, q, th))
#> [1] "mu_mu"       "sigma_sigma" "nu_nu"       "mu_sigma"    "mu_nu"      
#> [6] "sigma_nu"   
```
