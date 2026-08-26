# Is This Distribution Already Truncated?

Reports whether an object belongs to either truncated class.
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
asks before wrapping, and collapses a nested truncation into a single
object over the intersection of the two intervals. Nesting would be
correct but would pay the quadrature cost of
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
twice for a law one truncation already describes.

## Usage

``` r
is_truncated(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`. Any other input returns
  `FALSE`; nothing raises.

## Value

A single logical.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
which uses it, and
[TruncatedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedContinuousDistrib.md)
and
[TruncatedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedDiscreteDistrib.md),
the two classes it recognizes.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

c(truncated = distributions7:::is_truncated(tn),
  parent = distributions7:::is_truncated(gaussian1_distrib()))
#> truncated    parent 
#>      TRUE     FALSE 

# Nesting collapses: one object, over the intersection.
t2 <- truncated(truncated(gaussian1_distrib(), lower = -1), upper = 2)
c(lower = t2@lower, upper = t2@upper)
#> lower upper 
#>    -1     2 
distributions7:::is_truncated(t2@parent_distrib)
#> [1] FALSE
```
