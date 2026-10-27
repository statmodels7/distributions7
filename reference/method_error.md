# Error Message for a Method That Is Not an Object

Builds the message of
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for a `method` that is neither a
[FisherScoring](https://statmodels7.github.io/distributions7/reference/FisherScoring.md)
object nor an optimizer object.

## Usage

``` r
method_error(method)
```

## Arguments

- method:

  The rejected value.

## Value

A single character string.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
