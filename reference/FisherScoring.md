# S7 Class for Fisher Scoring Specifications

The class
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
returns and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
accepts through its `method` argument. It holds the two things a Fisher
scoring run needs that an ordinary optimizer cannot carry: how the
expected information is to be obtained when the family has no closed
form for it, and the stopping rule and iteration budget the run is to
use.

This page documents the raw S7 constructor, which validates nothing.
Build one with
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which checks all four arguments.

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

  How the expectation is approximated when the family has no closed-form
  expected information: one of `"bartlett"`, `"integrate"`, `"mc"` or
  `"opg"`. Ignored by a family that has one, and
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  rejects the call in that case.

- nsim:

  Monte Carlo sample size, a single positive number, read only when
  `approx = "mc"` is in force.

- criterion:

  A stopping rule from optimizers7, or `NULL` to take the default of
  [`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html).

- maxit:

  An iteration limit, a single positive number, or `NULL` for the same
  default.

## Value

An S7 object of class `FisherScoring` with those four properties.

## See also

[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which builds and validates one;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which consumes it.

## Examples

``` r
fs <- fisher_scoring(approx = "mc", nsim = 2000)
fs
#> Fisher scoring
#>   expected information: mc (nsim = 2000)  [ignored when the family has a closed form]
c(approx = fs@approx, nsim = fs@nsim)
#> approx   nsim 
#>   "mc" "2000" 

# The defaults leave the stopping rule and the budget to the fit.
d <- fisher_scoring()
c(criterion = is.null(d@criterion), maxit = is.null(d@maxit))
#> criterion     maxit 
#>      TRUE      TRUE 
```
