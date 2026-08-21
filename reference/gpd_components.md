# Derivative Components of the Generalized Pareto

The components of \\\partial^{a+b}\ell/\partial\sigma^a\partial\xi^b\\
at any order from one to four.

## Usage

``` r
gpd_components(y, theta, order, cut = 0.2, threads = 1L)
```

## Arguments

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- order:

  The derivative order, 1 to 4.

- cut:

  The value of \\\lvert\xi z\rvert\\ below which the series is used.

- threads:

  How many threads the series kernel may use.

## Value

A named list of component vectors, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

The shape direction goes through the series of \\\log(1+\xi z)/\xi\\
wherever \\\lvert\xi z\rvert\\ is below `cut`, which is the region where
the Leibniz form's terms of size \\\xi^{-(b+1)}\\ cancel against each
other; elsewhere the Leibniz form is used directly. The threshold is
exposed so that a test can force either route where both are accurate
and compare them.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md),
[`fdb2`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
