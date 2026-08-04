# Fisher Scoring as an Object

The S7 class of Fisher scoring specifications, returned by
[`fisher_scoring`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
and passed to
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
through its `method` argument.

## Usage

``` r
FisherScoring(
  approx = character(0),
  nsim = integer(0),
  criterion = NULL,
  maxit = NULL
)
```

## Arguments

- approx:

  How the expectation is approximated when the distribution has no
  closed-form expected information.

- nsim:

  Monte Carlo sample size, used when `approx = "mc"`.

- criterion:

  The stopping rule, or `NULL` to take the fit's.

- maxit:

  The iteration limit, or `NULL` to take the fit's.

## Value

An object of class `FisherScoring`.

## See also

[`fisher_scoring`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)

## Examples

``` r
fs <- fisher_scoring(approx = "mc", nsim = 2000)
S7::S7_inherits(fs, FisherScoring)
#> [1] TRUE
fs@approx
#> [1] "mc"
```
