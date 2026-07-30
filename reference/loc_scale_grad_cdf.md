# Location-Scale CDF Gradient

The
[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
body shared by the location-scale families:
[`loc_scale_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
at order 1, put on the requested tail and scale.

## Usage

``` r
loc_scale_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE)
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

A named list of gradient component vectors.
