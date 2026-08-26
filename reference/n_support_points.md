# Number of Points in a Discrete Support

Returns how many points the distribution's support contains, read off
its `bounds`, and `Inf` when the support is unbounded above. It is what
[`check_support_is_rich_enough()`](https://statmodels7.github.io/distributions7/reference/check_support_is_rich_enough.md)
counts against the parameters a zero wrapper would spend.

## Usage

``` r
n_support_points(distrib)
```

## Arguments

- distrib:

  A `distrib` object, whose `bounds` are read.

## Value

A single number: the count of integer points between the bounds
inclusive, or `Inf`.

## See also

[`check_support_is_rich_enough()`](https://statmodels7.github.io/distributions7/reference/check_support_is_rich_enough.md),
the consumer, and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
for the counting rule.

## Examples

``` r
# A count family has an unbounded support.
distributions7:::n_support_points(poisson_distrib())
#> [1] Inf

# A binomial has size + 1 points, and a Bernoulli has two.
c(binomial5 = distributions7:::n_support_points(binomial_distrib(size = 5)),
  bernoulli = distributions7:::n_support_points(bernoulli_distrib()))
#> binomial5 bernoulli 
#>         6         2 
```
