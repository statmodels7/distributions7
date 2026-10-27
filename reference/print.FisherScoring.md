# Print a Fisher Scoring Specification

Writes the method that obtains the expected information, with the Monte
Carlo sample size when that is the method, and then any setting that
overrides the fit's: a stopping rule set here is reported as set, and an
iteration limit by its value. For `method = NULL` the line states that
the closed form of the family is used where there is one.

## Arguments

- x:

  A
  [`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md)
  object.

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. The function is called for the output that it writes.

## See also

[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for what each line reports.

## Examples

``` r
fisher_scoring()
#> Fisher scoring
#>   expected information: closed form where the family has one, otherwise the outer product of the scores
fisher_scoring(method = expect_mc(nsim = 2000))
#> Fisher scoring
#>   expected information: by Monte Carlo (nsim = 2000)
fisher_scoring(criterion = optimizers7::crit_grad(1e-9), maxit = 50)
#> Fisher scoring
#>   expected information: closed form where the family has one, otherwise the outer product of the scores
#>   criterion: set here
#>   maxit: 50
```
