# Print an Expectation Method

Writes one line naming the method and, for
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md),
the number of draws, and for
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
the identity when it is not the default.

## Arguments

- x:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. The function is called for the output that it writes.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what each method computes.

## Examples

``` r
expect_opg()
#> Expectation by the outer product of the scores
expect_mc(nsim = 500)
#> Expectation by Monte Carlo (nsim = 500)
expect_integrate(identity = "score")
#> Expectation by integration against the density (identity = score)
```
