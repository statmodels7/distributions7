# The Rows of a Family's Constants Inside an Integrand

Prepares the integrand of
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
and
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
for a family whose constants vary by observation (a binomial's `size`).
The observation index travels to the integrand as the further argument
`const_row` of
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
which indexes it like a parameter, and `at()` returns the family
restricted to those rows with
[`distrib_at_rows()`](https://statmodels7.github.io/distributions7/reference/distrib_at_rows.md).
For a family without such constants `dots` is empty and `at()` returns
the family unchanged.

## Usage

``` r
expectation_rows(distrib)
```

## Arguments

- distrib:

  A univariate family.

## Value

A list with `dots`, the further arguments for
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
and `at`, a function of those arguments that returns the family.

## Examples

``` r
r <- distributions7:::expectation_rows(binomial_distrib(size = c(5, 10)))
r$at(const_row = c(2, 2, 1))@size
#> [1] 10 10  5
```
