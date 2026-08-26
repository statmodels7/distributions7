# Does This Distribution Already Model a Probability of Zero

Answers whether the object is already one of the two zero wrappers,
which is how each of them refuses to be applied on top of the other.
Stacking them leaves only their combination identified, so the check is
what turns a model with a flat direction into an error at construction.

## Usage

``` r
is_zero_wrapper(distrib)
```

## Arguments

- distrib:

  A `distrib` object.

## Value

`TRUE` for a zero-inflated or a zero-adjusted object of either kind,
`FALSE` otherwise.

## See also

[`check_not_stacked()`](https://statmodels7.github.io/distributions7/reference/check_not_stacked.md),
which reports the refusal, and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for what stacking would cost.

## Examples

``` r
c(plain = distributions7:::is_zero_wrapper(poisson_distrib()),
  inflated = distributions7:::is_zero_wrapper(
    zero_inflated(poisson_distrib())),
  adjusted = distributions7:::is_zero_wrapper(
    zero_adjusted(poisson_distrib())))
#>    plain inflated adjusted 
#>    FALSE     TRUE     TRUE 
```
