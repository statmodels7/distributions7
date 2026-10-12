# The Log-Probabilities on the Additive Log-Ratio Chart

Computes \\\log p_j = \eta_j - \mathrm{LSE}(\eta_1, \ldots, \eta\_{K-1},
0)\\, with \\\eta_K = 0\\, the log-sum-exp shifted by its largest
argument. A probability that underflows keeps its logarithm.

## Usage

``` r
mn_log_probs(s, eta)
```

## Arguments

- s:

  A
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  parameter.

- eta:

  Its free vector.

## Value

A numeric vector of length `s@n_free + 1`.

## See also

[`mn_is_alr()`](https://statmodels7.github.io/distributions7/reference/mn_is_alr.md)
and
[`mn_cumulant()`](https://statmodels7.github.io/distributions7/reference/mn_cumulant.md).

## Examples

``` r
distributions7:::mn_log_probs(parameters7::simplex(3), c(-800, 0))
#> [1] -800.6931472   -0.6931472   -0.6931472
```
