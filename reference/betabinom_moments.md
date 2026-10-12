# Moments of the Beta-Binomial in Closed Form

Returns the mean, the variance, the skewness or the excess kurtosis of
the beta-binomial with \\n\\ trials, mean proportion \\\mu\\ and
dispersion \\\sigma = 1/(\alpha + \beta)\\, each in a closed form that
forms no difference of large raw moments. With \\g = \mu(1-\mu)\\,
\$\$E\[Y\] = n\mu,\qquad \mathrm{Var}(Y) = ng\\\frac{1 + n\sigma}{1 +
\sigma},\$\$ \$\$\gamma_1 = (1 - 2\mu)\\\frac{1 + 2n\sigma}{1 + 2\sigma}
\sqrt{\frac{1 + \sigma}{gn(1 + n\sigma)}},\$\$ \$\$\gamma_2 =
\frac{6n\sigma(1 + n\sigma)\\\[(1 - 6g)(1 + \sigma) + g\] + (1 +
\sigma)(1 - 6g - \sigma)}{gn(1 + 2\sigma)(1 + 3\sigma)(1 +
n\sigma)}.\$\$ At \\\sigma = 0\\ these are the binomial's, \\\gamma_2\\
reducing to \\(1 - 6g)/(gn)\\. The excess kurtosis is written as one
quotient, so that it is not formed as a kurtosis near 3 minus 3, which
loses the digits of a value of order \\1/n\\.

## Usage

``` r
betabinom_moments(mu, sigma, n, what)
```

## Arguments

- mu:

  The mean proportion, a numeric vector in \\(0, 1)\\.

- sigma:

  The dispersion, a positive numeric vector.

- n:

  The number of trials, a single positive whole number.

- what:

  One of `"mean"`, `"variance"`, `"skewness"` or `"kurtosis"`.

## Value

A numeric vector of the recycled length of `mu` and `sigma`.

## See also

[`mean.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom1Distrib.md),
[`variance.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom1Distrib.md),
[`skewness.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom1Distrib.md)
and
[`kurtosis.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md),
and the methods of
[BetaBinom2Distrib](https://statmodels7.github.io/distributions7/reference/BetaBinom2Distrib.md),
which call it.

## Examples

``` r
# At a small dispersion the excess kurtosis tends to the binomial's.
distributions7:::betabinom_moments(0.3, 1e-12, 1e6, "kurtosis")
#> [1] -1.238095e-06
(1 - 6 * 0.21) / (0.21 * 1e6)
#> [1] -1.238095e-06
```
