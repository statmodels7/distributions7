# Moments of a Reparametrized Distribution

The parent's, at the mapped parameters. A reparametrization does not
change the law, so it does not change a moment: what changes is the
coordinates the moment is computed from. Delegating keeps the parent's
closed forms instead of falling through to a quadrature.

## Usage

``` r
reparam_mean(x, theta, ...)

reparam_variance(x, theta, ...)

reparam_skewness(x, theta, ...)

reparam_kurtosis(x, theta, ...)
```

## Arguments

- x:

  A reparametrized distribution.

- theta:

  A named list of the new parameters.

- ...:

  Passed to the parent.

## Value

A numeric vector.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
