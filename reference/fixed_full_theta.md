# Splice the Fixed Values Back Into a Full Parameter List

Combines the wrapper's free `theta` with its fixed values into the full
parameter list the PARENT expects, in the parent's own order. Every
method of every `Fixed*` class begins with this call and then delegates,
so it is the single point at which the two parameter sets are
reconciled.

## Usage

``` r
fixed_full_theta(distrib, theta)
```

## Arguments

- distrib:

  A `FixedContinuousDistrib`, `FixedDiscreteDistrib` or
  `FixedMultivariateDistrib` object.

- theta:

  A named list of the FREE parameters, already aligned. A fully-fixed
  wrapper takes [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list of the parent's parameters, complete and in the parent's
order.

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
for the wrapper and
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
for what the methods do with the result.

## Examples

``` r
d <- fixed(gaussian1_distrib(), mu = 0)
str(distributions7:::fixed_full_theta(d, list(sigma = 2)))
#> List of 2
#>  $ mu   : num 0
#>  $ sigma: num 2

# Which is exactly what the parent is then called at.
all.equal(distrib_pdf(d, c(-1, 1), list(sigma = 2)),
          distrib_pdf(gaussian1_distrib(), c(-1, 1),
                      distributions7:::fixed_full_theta(d, list(sigma = 2))))
#> [1] TRUE
```
