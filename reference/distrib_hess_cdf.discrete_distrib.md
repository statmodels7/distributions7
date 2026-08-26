# Log-CDF Hessian for Discrete Distributions

Exact, by the same finite sum as the gradient: at second order the
summand is \\\ell^{(ij)} + \ell^{(i)}\ell^{(j)}\\, the Bartlett lemma's
expansion of \\\partial^2 f/f\\. The first-order part comes from
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
so a family's own simplification is used where it has one.

## Arguments

- distrib:

  A `discrete_distrib` object.

- q:

  A numeric vector of quantiles. Values below the support give zero.

- theta:

  A named list of parameters on the parameter scale.

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

## See also

[`distrib_grad_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md)
for the first order;
[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
for the sum;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).

## Examples

``` r
# A negative binomial: two parameters, so three Hessian components.
d <- negbin2_distrib()
distrib_hess_cdf(d, c(2, 5), list(mu = 3, theta = 2))
#> $mu_mu
#> [1]  0.01841761 -0.02370989
#> 
#> $theta_theta
#> [1]  0.041640393 -0.007702265
#> 
#> $mu_theta
#> [1] -0.060820940 -0.007018114
#> 
```
