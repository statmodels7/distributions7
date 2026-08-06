# Derivative Components of the Beta-Binomial in Mean and Dispersion

The components of any order from one to four, obtained by carrying the
shape parametrization's derivatives through the map \\(a, b) =
(\mu/\sigma, (1-\mu)/\sigma)\\.

## Usage

``` r
betabinom1_components(distrib, y, theta, order)
```

## Arguments

- distrib:

  A
  [`BetaBinom1Distrib`](https://statmodels7.github.io/distributions7/reference/BetaBinom1Distrib.md)
  object.

- y:

  A numeric vector of counts.

- theta:

  A list containing `mu` and `sigma`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of component vectors, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
[`chain_derivatives`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
