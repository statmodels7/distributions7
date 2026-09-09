# One Skew t Derivative in the Location, Scale and Shape

Reads \\\partial^a\_\mu \partial^b\_\sigma \partial^c\_\alpha \ell\\ off
the table
[`skewt_msa_tower()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_tower.md)
builds, by the \\P\_{a,b}\\ sum that function's page derives.

## Usage

``` r
skewt_msa_component(tw, sigma, a, b, c)
```

## Arguments

- tw:

  A list from
  [`skewt_msa_tower()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_tower.md).

- sigma:

  The scale, of length 1 or of the length of the response.

- a, b, c:

  Non-negative integers, the number of times the derivative is taken in
  \\\mu\\, in \\\sigma\\ and in \\\alpha\\. Their sum is the order.

## Value

A numeric vector, one value per observation.

## Details

A term whose coefficient \\\prod\_{k=i}^{b-1}(a+k)\\ vanishes is dropped
rather than multiplied, because the factor it would multiply is the
entry `"0_0"`, which the table does not hold: in R `0 * NA` is `NA`.

## See also

[`skewt_msa_tower()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_tower.md)
for the table and the derivation, and
[`skewt_msa_derivs()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_derivs.md),
which loops this over a whole order.

## Examples

``` r
tw <- distributions7:::skewt_msa_tower(c(-0.4, 1.2), 0, 1, 0.7, 8)
distributions7:::skewt_msa_component(tw, 1, a = 2L, b = 1L, c = 0L)
#> [1] 2.257940 0.811152
```
