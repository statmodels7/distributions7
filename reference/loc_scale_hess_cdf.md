# Location-Scale CDF Hessian

The
[`distrib_hess_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
body shared by the location-scale families.

## Usage

``` r
loc_scale_hess_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE)
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
