# Pseudo-Huber Random Generation, Standard-Deviation Parametrization

Draws `n` variates as the normal variance mixture of
[`distrib_rng.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuberDistrib.md),
at the scale \\\sigma_1 = \sigma / \sqrt{R(\nu)}\\ of that family.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- n:

  The number of draws.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- ...:

  Unused.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_rng.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuberDistrib.md).

## Examples

``` r
set.seed(1)
x <- distrib_rng(pseudohuber2_distrib(), 200,
                 list(mu = 0, sigma = 2, nu = 3))
sd(x)
#> [1] 2.121772
```
