# Reject the Composition of Two Zero Wrappers

Signals an error when the parent already models the probability of a
zero, which is how both wrappers refuse to stack. The models this
rejects are well defined and inestimable, so the constructor is the only
place the problem can be caught: nothing at run time reports it.

## Usage

``` r
check_not_stacked(distrib, fun, param)
```

## Arguments

- distrib:

  The candidate parent, a `distrib` object.

- fun:

  The name of the calling constructor, a single string, quoted in the
  message.

- param:

  The name of the parameter it would add, a single string. Not currently
  placed in the message; it is carried for the caller's own record.

## Value

`NULL`, invisibly, when the parent is not a zero wrapper. Otherwise it
signals an error naming the parent.

## Details

Two zero parameters cannot both be identified. Zero-adjusting a
zero-inflated parent truncates the zero away, which cancels the
\\(1-\zeta)\\ factor between the numerator and the normalizing constant,
so \\\zeta\\ leaves the likelihood entirely and its score is IDENTICALLY
zero. The other order keeps only the combination \\\zeta +
(1-\zeta)\pi\\. Either way an optimizer wanders along a flat ridge, the
mass function sums to one throughout, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
passes.

## See also

[`is_zero_wrapper()`](https://statmodels7.github.io/distributions7/reference/is_zero_wrapper.md)
for the test, and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
the two callers.

## Examples

``` r
# A plain parent passes silently.
distributions7:::check_not_stacked(poisson_distrib(),
                                   "zero_inflated()", "zi")

# A wrapped one does not, and the two orders are refused alike.
try(zero_inflated(zero_adjusted(poisson_distrib())))
#> Error : zero_inflated() cannot wrap 'zero-adjusted poisson', which already models the probability of a zero.
#>   Stacking the two leaves only their combination identified: the second
#>   parameter has an identically zero score, and any optimizer will wander
#>   along that ridge. Apply exactly one zero wrapper to a plain distribution.
try(zero_adjusted(zero_inflated(poisson_distrib())))
#> Error : zero_adjusted() cannot wrap 'zero-inflated poisson', which already models the probability of a zero.
#>   Stacking the two leaves only their combination identified: the second
#>   parameter has an identically zero score, and any optimizer will wander
#>   along that ridge. Apply exactly one zero wrapper to a plain distribution.
```
