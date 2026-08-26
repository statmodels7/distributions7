# Location-Scale Third and Fourth Log-CDF Derivatives

Builds the
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
or
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
body that the location-scale families register:
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
at every order up to the one wanted, put on the requested tail and scale
by
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).
Five families use it, the Gaussian, the logistic, the Cauchy and the
Laplace here and the Gumbel in `cdf_mapped_higher.R`.

## Usage

``` r
loc_scale_deriv_cdf_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function of `(distrib, q, theta, lower.tail, log, ...)` suitable for
registering as an S7 method on either generic, returning a named list of
numeric vectors of that order.

## Details

Every lower order is computed because the log-scale conversion reads
them all. The function is a factory so that the order can be closed
over, `force(order)` being what keeps the two registrations from sharing
one value.

## See also

[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
for the formulas;
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
for the generic;
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
and
[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
for the orders below.
