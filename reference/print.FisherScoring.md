# Print a Fisher Scoring Specification

Writes the strategy for the expected information, with the Monte Carlo
sample size when that is the strategy, and then any setting that
overrides the fit's: a stopping rule set here is reported as set, and an
iteration limit by its value. The line about the strategy carries the
reminder that a family with a closed-form expectation ignores it.

## Arguments

- x:

  A
  [`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md)
  object.

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. Called for the output it writes.

## See also

[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for what each line reports.

## Examples

``` r
fisher_scoring()
#> Fisher scoring
#>   expected information: opg  [ignored when the family has a closed form]
fisher_scoring(approx = "mc", nsim = 2000)
#> Fisher scoring
#>   expected information: mc (nsim = 2000)  [ignored when the family has a closed form]
fisher_scoring(criterion = optimizers7::crit_grad(1e-9), maxit = 50)
#> Fisher scoring
#>   expected information: opg  [ignored when the family has a closed form]
#>   criterion: set here
#>   maxit: 50
```
