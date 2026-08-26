# The Dispersion a Poisson-Inverse Gaussian Alpha Implies

Converts the orthogonal parametrization's \\\alpha\\ into the dispersion
\\\sigma\\ of
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md):
\$\$\sigma = \dfrac{\mu + \sqrt{\mu^2 + \alpha^2}}{\alpha^2},\$\$ the
positive root of \\\alpha^2\sigma^2 - 2\mu\sigma - 1 = 0\\.

## Usage

``` r
pig2_sigma(mu, alpha)
```

## Arguments

- mu:

  A numeric vector of means, strictly positive. Nothing is validated.

- alpha:

  A numeric vector of Bessel arguments, strictly positive.

## Value

A numeric vector of dispersions, of the length of the recycled inputs,
strictly positive.

## Details

The relation is the inverse of \\\alpha = \sqrt{1 +
2\sigma\mu}/\sigma\\. Written this way there is no cancellation anywhere
in the domain: both \\\mu\\ and \\\sqrt{\mu^2 + \alpha^2}\\ are
positive, so the numerator is a sum of positive terms. The other root of
the quadratic is negative and is discarded.

Only
[`distrib_rng.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig2Distrib.md)
calls this. The derivative methods do not: the compiled kernel takes
\\\alpha\\ as a variable of its own, so the map is never differentiated.

## See also

[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the parametrization,
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the one it maps onto, and
[`distrib_rng.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig2Distrib.md)
for its only caller.

## Examples

``` r
# The map, and the round trip back through alpha = sqrt(1 + 2 sigma mu)/sigma.
sg <- distributions7:::pig2_sigma(3, 3.010399)
c(sigma = sg, alpha_back = sqrt(1 + 2 * sg * 3) / sg)
#>      sigma alpha_back 
#>  0.7999998  3.0103990 

# A large alpha is a small dispersion: the family tends to the Poisson.
distributions7:::pig2_sigma(3, c(1, 3, 30, 300))
#> [1] 6.162277660 0.804737854 0.036832919 0.003366833
```
