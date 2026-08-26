# Elastic-Net Random Generation

Draws by the inverse transform,
`distrib_quantile(distrib, runif(n), theta)`. The quantile function of
this family is closed form, so the transform is exact and costs one
uniform per draw; the base class's ratio-of-uniforms fallback is not
needed.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of length `n`; a component of length 1
  is recycled.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws, symmetric about `mu`.

## See also

[`distrib_quantile.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.EnetDistrib.md),
which it inverts through,
[`variance.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.EnetDistrib.md)
for the moment the draws reproduce, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 0, lambda = 2, alpha = 0.5)

set.seed(71)
x <- distrib_rng(d, 2e5, th)

# Two moments against their closed forms, and the symmetry.
rbind(sample = c(mean(x), var(x), mean((x - mean(x))^3) / sd(x)^3),
      theory = c(mean(d, th), variance(d, th), skewness(d, th)))
#>                 [,1]      [,2]         [,3]
#> sample -0.0005352337 0.4778600 -0.008027417
#> theory  0.0000000000 0.4748647  0.000000000

# As alpha approaches one the draws become Laplace.
set.seed(72)
xl <- distrib_rng(d, 1e5, list(mu = 0, lambda = 2, alpha = 1 - 1e-10))
c(sample_var = var(xl), laplace_var = 2 / 2^2)
#>  sample_var laplace_var 
#>   0.4889258   0.5000000 
```
