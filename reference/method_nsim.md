# Monte Carlo Sample Size of an Expectation Method

Returns the number of draws carried by an
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
object, and 10000, the default of
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md),
for any other method. The multivariate third and fourth derivatives take
their expectation by sampling whatever the class of `method` and read
the sample size here.

## Usage

``` r
method_nsim(method)
```

## Arguments

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  or `NULL`.

## Value

A single positive number.

## See also

[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md),
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
