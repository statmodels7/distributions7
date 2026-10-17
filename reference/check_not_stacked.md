# Reject the Composition of Two Zero Wrappers

Signals an error when the parent already models the probability of a
zero, which is how both constructors reject a parent that is already a
zero wrapper. The models that this function rejects are well defined and
inestimable, and the problem is not reported at run time, so the
constructor is the only place where it can be caught.

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

  The name of the parameter that the calling constructor would add, a
  single string. An error is signaled when the parent already has a
  parameter with this name.

## Value

`NULL`, invisibly, when the parent is not a zero wrapper and has no
parameter named `param`. If the parent is a zero wrapper, the function
signals an error that names the parent; if the parent already has a
parameter named `param`, it signals an error that names that parameter.

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
"zero_inflated", "zi")

# A wrapped one does not, and the two orders are rejected alike.
try(zero_inflated(zero_adjusted(poisson_distrib())))
#> Error : zero_inflated() cannot wrap 'zero-adjusted poisson', which already models the probability of a zero.
#>   Stacking the two leaves only a combination of their parameters identified,
#>   so the likelihood is flat along a direction and any optimizer will wander
#>   along that ridge. Apply exactly one zero wrapper to a plain distribution.
try(zero_adjusted(zero_inflated(poisson_distrib())))
#> Error : zero_adjusted() cannot wrap 'zero-inflated poisson', which already models the probability of a zero.
#>   Stacking the two leaves only a combination of their parameters identified,
#>   so the likelihood is flat along a direction and any optimizer will wander
#>   along that ridge. Apply exactly one zero wrapper to a plain distribution.
```
