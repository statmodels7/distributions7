# Reject a Model With More Parameters Than the Support Can Distinguish

Signals an error when the parent's support has too few points for the
zero wrapper's extra parameter to be identified. A discrete distribution
on \\k\\ points has \\k - 1\\ free probabilities, and a wrapper spends
`n_params + 1` of them, so \\k \ge n\_{params} + 2\\ is necessary. The
rule is the SAME for both wrappers.

## Usage

``` r
check_support_is_rich_enough(distrib, fun)
```

## Arguments

- distrib:

  The candidate parent, a `distrib` object.

- fun:

  The name of the calling constructor, a single string, quoted in the
  message.

## Value

`NULL`, invisibly, when the support is large enough. Otherwise it
signals an error giving the point count, the parameter count and the
count required.

## Details

What it rules out is exactly the Bernoulli and
`binomial_distrib(size = 1)`. Zero-adjusting a Bernoulli leaves the
truncated part on the single point \\\\1\\\\ and `mu` disappears: the
mass function is literally the same at `mu = 0.2` and at `mu = 0.9`.
None of this is visible at run time, the mass summing to one and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
passing, so the constructor is the only place it can be caught.

A large support is NECESSARY without being sufficient. With a mean large
enough that the parent puts almost no mass at zero, the extra parameter
is weakly identified whatever the support size, which is a question
about the data and one this check cannot ask.

## See also

[`n_support_points()`](https://statmodels7.github.io/distributions7/reference/n_support_points.md)
for the count, and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
the two callers.

## Examples

``` r
# A count family has room.
distributions7:::check_support_is_rich_enough(poisson_distrib(),
                                              "zero_inflated()")

# A Bernoulli does not, and neither does a one-trial binomial.
try(zero_inflated(bernoulli_distrib()))
#> Error : zero_inflated() cannot wrap 'bernoulli': its support has 2 points, so the family has 1 free
#>   probability, while the wrapped distribution would have 2 parameters. They are not
#>   identified -- different parameter values give exactly the same distribution.
#>   A support of at least 3 points is required.
try(zero_inflated(binomial_distrib(size = 1)))
#> Error : zero_inflated() cannot wrap 'binomial': its support has 2 points, so the family has 1 free
#>   probability, while the wrapped distribution would have 2 parameters. They are not
#>   identified -- different parameter values give exactly the same distribution.
#>   A support of at least 3 points is required.

# Five trials is enough: six points against two parameters.
zero_inflated(binomial_distrib(size = 5))@params
#> [1] "mu" "zi"
```
