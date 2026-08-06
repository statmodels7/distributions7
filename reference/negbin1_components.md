# Derivative Components of NB1

The components of \\\partial^{a+b}\ell/\partial\mu^a\partial\theta^b\\
at any order from one to four, from the sparse form of the
log-likelihood in the size \\r = \mu/\theta\\.

## Usage

``` r
negbin1_components(y, theta, order)
```

## Arguments

- y:

  A numeric vector of counts.

- theta:

  A list containing `mu` and `theta`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of component vectors, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
