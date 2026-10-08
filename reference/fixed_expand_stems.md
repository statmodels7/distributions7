# Expand the Stem of a Vector Parameter in fixed()

Replaces a name that is not a parameter but is the stem of the
parameters `stem1, ..., stemk` by those parameters, recycling a single
value or taking one value per coordinate. It is what lets
`fixed(mvgaussian1_distrib(3), mu = 0)` hold the whole mean at zero.

## Usage

``` r
fixed_expand_stems(fix, params)
```

## Arguments

- fix:

  The named list of values passed to
  [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md).

- params:

  The parameter names of the distribution being fixed.

## Value

The list with every stem expanded into its coordinates. A name that is a
parameter, or matches no stem, is returned unchanged for
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
to check.

## Examples

``` r
distributions7:::fixed_expand_stems(list(mu = 0, sigma = 1),
                                    c("mu1", "mu2", "sigma"))
#> $mu1
#> [1] 0
#> 
#> $mu2
#> [1] 0
#> 
#> $sigma
#> [1] 1
#> 
distributions7:::fixed_expand_stems(list(mu = c(1, 2)), c("mu1", "mu2"))
#> $mu1
#> [1] 1
#> 
#> $mu2
#> [1] 2
#> 
```
