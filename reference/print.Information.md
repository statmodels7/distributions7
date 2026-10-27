# Print an Information Object

Writes one line naming the matrix and, for
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md),
the method that obtains the expectation.

## Arguments

- x:

  An object inheriting from
  [Information](https://statmodels7.github.io/distributions7/reference/Information.md).

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. The function is called for the output that it writes.

## See also

[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)

## Examples

``` r
observed_info()
#> Observed information
expected_info(method = expect_mc(nsim = 500))
#> Expected information (expectation by Monte Carlo (nsim = 500))
```
