# The Beta-Binomial's Expected Higher Derivatives in Mean and Dispersion

Returns the expectation of every derivative of one order of the
beta-binomial log-mass in \\(\mu, \sigma)\\ as an exact finite sum over
the support \\\\0, \dots, n\\\\, the mass times the derivative at each
count, through the same summation as
[`betabinom2_expected()`](https://statmodels7.github.io/distributions7/reference/betabinom2_expected.md).
Parameters that vary by observation give one sum per observation.

## Usage

``` r
betabinom1_expected(distrib, y, theta, order)
```

## Arguments

- distrib:

  A
  [`BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/BetaBinom1Distrib.md)
  object.

- y:

  A numeric vector; only its length is read.

- theta:

  A named list with components `mu` and `sigma`.

- order:

  The derivative order, 3 or 4.

## Value

A named list of numeric vectors keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them, each of the recycled length of `y` and the parameters.

## See also

[`betabinom1_components()`](https://statmodels7.github.io/distributions7/reference/betabinom1_components.md)
for the summand,
[`distrib_deriv3.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md)
and
[`distrib_deriv4.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom1Distrib.md),
which call this.
