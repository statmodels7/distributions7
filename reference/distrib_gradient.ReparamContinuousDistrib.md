# Gradient of a Reparametrized Distribution

The parent's score carried by the Jacobian of the map, \\\ell^{(a)} =
\sum_i \ell^{(i)} \partial\theta_i/\partial\psi_a\\.

## Usage

``` r
reparam_gradient(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response.

- theta:

  A named list of the new parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of first derivatives.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
