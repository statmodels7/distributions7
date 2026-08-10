# The Mixed Grid of the Lognormal

The gaussian's components at \\t = \log y\\, carried by the Jacobian of
the transformation.

## Usage

``` r
lognormal_theta_chain(y, theta, order, second)
```

## Arguments

- y:

  A numeric vector of observations.

- theta:

  A named list containing `mu` and `sigma2`.

- order:

  `1` for the block of
  [`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` for that of
  [`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).

- second:

  Whether the second-order theta components are wanted.

## Value

A named list, keyed by parameter or by parameter pair.

## See also

[`distrib_grad_y_hess`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
