# Whether a Simplex Parameter Is the Additive Log-Ratio Chart

Returns `TRUE` when the probability parameter is a
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html),
whose free values are the additive log-ratios against the last category.
On that chart the log-probabilities are \\\eta_j - \mathrm{LSE}(\eta)\\
and their derivatives are cumulants of the categorical law, which
[`mn_log_probs()`](https://statmodels7.github.io/distributions7/reference/mn_log_probs.md)
and
[`mn_cumulant()`](https://statmodels7.github.io/distributions7/reference/mn_cumulant.md)
compute without dividing by a probability.

## Usage

``` r
mn_is_alr(s)
```

## Arguments

- s:

  The probability parameter of a
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
  or a
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

## Value

A single logical.

## See also

[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md)
and
[`distrib_hessian.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MultinomialDistrib.md).

## Examples

``` r
distributions7:::mn_is_alr(parameters7::simplex(3))
#> [1] TRUE
```
