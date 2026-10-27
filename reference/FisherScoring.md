# S7 Class for Fisher Scoring Specifications

This class is the one that
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
returns and that
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
accepts through its `method` argument. It holds the two things that a
Fisher scoring run needs and that an ordinary optimizer cannot carry:
the method that obtains the expected information, `NULL` or an object of
class
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
and the stopping rule and iteration budget that the run uses.

This page documents the raw S7 constructor, which validates nothing.
Build one with
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which checks all four arguments.

## Usage

``` r
FisherScoring(method = NULL, criterion = NULL, maxit = NULL)
```

## Arguments

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
  that selects how the expected information is obtained. `NULL` takes
  the closed form of the family where it has one and
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  otherwise; an object is used for every family, also for one that has
  the closed form.

- criterion:

  A stopping rule from optimizers7, or `NULL` to take
  [`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)
  at its default tolerance.

- maxit:

  An iteration limit, a single positive number, or `NULL` for the
  default of
  [`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html).

## Value

An S7 object of class `FisherScoring` with the properties `method`
(`NULL` or an object inheriting from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)),
`criterion` and `maxit`.

## See also

[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which builds and validates one;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which consumes it;
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
for the class of the `method` property.

## Examples

``` r
fs <- fisher_scoring(method = expect_mc(nsim = 2000))
fs
#> Fisher scoring
#>   expected information: by Monte Carlo (nsim = 2000)
fs@method@nsim
#> [1] 2000
is.null(fisher_scoring()@method)
#> [1] TRUE

# The defaults leave the stopping rule and the budget to the fit.
d <- fisher_scoring()
c(criterion = is.null(d@criterion), maxit = is.null(d@maxit))
#> criterion     maxit 
#>      TRUE      TRUE 
```
