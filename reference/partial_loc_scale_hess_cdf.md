# CDF Hessian When Only Some Parameters Are Location-Scale

The second-order companion of
[`partial_loc_scale_grad_cdf`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md):
the three components in the location and scale in closed form, the
components involving a shape parameter by finite differences.

## Usage

``` r
partial_loc_scale_hess_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

## Value

A named list of Hessian component vectors.

## See also

[`loc_scale_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
