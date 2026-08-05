# Observed Hessian of a Reparametrized Distribution

The second-order chain rule, which keeps the term in the parent's score
and the second derivative of the map: \\\ell^{(ab)} =
\sum\_{ij}\ell^{(ij)}h^i_a h^j_b + \sum_i \ell^{(i)} h^i\_{ab}\\.

## Usage

``` r
reparam_hessian(distrib, y, theta, scale = c("parameter", "link"), ...)
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

A named list of second derivatives.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
