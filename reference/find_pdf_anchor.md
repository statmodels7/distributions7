# Locate an Interior High-Density Point

Finds an approximate mode, used to split integrals, to scale
root-finding brackets and to recentre the ratio-of-uniforms kernel.

## Usage

``` r
find_pdf_anchor(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameters.

## Value

A single number, an interior point of high density.

## Details

The support is compactified to \\(0, 1)\\ and the maximum found by
repeatedly evaluating the log-density on a grid and keeping the two
cells around the largest value: for a unimodal density that bracket
provably still contains the mode, and it shrinks by a factor of 64 per
pass.

A golden-section search on the compactified scale is *not* accurate
enough here, which is worth knowing before anyone replaces this with
one. Its tolerance is expressed in the compactified variable, and the
derivative of the compactification can be enormous – under the tangent
map \\dy/dt\\ is of order \\y^2\\, so a default tolerance of about
`1e-4` became an error of 125 standard deviations for a density centred
at 1000. The grid refinement instead stops on the width of the bracket
measured in \\y\\.

## See also

[`find_lp_anchor`](https://statmodels7.github.io/distributions7/reference/find_lp_anchor.md),
[`rng_grou`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
