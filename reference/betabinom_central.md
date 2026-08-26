# Mean and Central Moments of a Beta-Binomial

The mean and the second, third and fourth central moments of a
beta-binomial, one value per parameter setting. It maps the
mean-dispersion parametrization onto the two shapes, \\a = \mu/\sigma\\
and \\b = (1-\mu)/\sigma\\, then chains
[`betabinom_factorial_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_factorial_moments.md)
into
[`central_from_factorial()`](https://statmodels7.github.io/distributions7/reference/central_from_factorial.md)
once per setting and stacks the results. The four moment methods of the
family share it.

## Usage

``` r
betabinom_central(mu, sigma, n)
```

## Arguments

- mu:

  The success probability, a numeric vector strictly inside \\(0,1)\\.

- sigma:

  The dispersion, a positive numeric vector. As it goes to zero the
  mixing beta concentrates and the family tends to a binomial.

- n:

  The number of trials, a non-negative whole number.

## Value

A named list with `mean` and the central moments `c2`, `c3` and `c4`,
each a numeric vector as long as the longer of `mu` and `sigma`.

## Notation

\\\mu \in (0,1)\\ is the success probability, \\\sigma \> 0\\ the
dispersion and \\n\\ the number of trials.

## See also

[`variance.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom1Distrib.md),
[`skewness.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom1Distrib.md)
and
[`kurtosis.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md)
for the consumers.

## Examples

``` r
# The mean is n mu, whatever the dispersion.
distributions7:::betabinom_central(0.3, 0.5, 10)$mean
#> [1] 3
```
