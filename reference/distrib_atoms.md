# Atoms of a Distribution

Returns the locations and probabilities of the point masses a
distribution places on individual values — the discrete part of a
*mixed* distribution, one that is neither purely continuous nor purely
discrete.
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
applied to a continuous distribution builds exactly such an object: a
point mass at zero next to a density.

## Usage

``` r
distrib_atoms(distrib, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- theta:

  A named list (or named numeric vector) of distribution parameters,
  with scalar entries.

- ...:

  Additional arguments passed to the specific method.

## Value

A list with components `y` (the locations) and `p` (their
probabilities), both numeric vectors of the same length, possibly of
length zero.

## Details

The default returns no atoms, which is the right answer for every
ordinary distribution: a continuous one has none, and a discrete one is
made of nothing else, so listing them would be pointless. The generic
exists for the case in between, where a routine written for densities
has to be told that part of the mass is not in the integral —
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
uses it to know that the density is expected to integrate to \\1 - \sum
p\\ rather than 1, and to keep its finite differences away from the
jumps.

## See also

[`distrib_pdf`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
[`distrib_quantile`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
[`distrib_rng`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)

## Examples

``` r
distrib_atoms(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#> $y
#> numeric(0)
#> 
#> $p
#> numeric(0)
#> 
distrib_atoms(zero_adjusted(gamma2_distrib()), list(mu = 2, sigma2 = 1, za = 0.3))
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 
```
