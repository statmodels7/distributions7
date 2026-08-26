# Locate an Interior High-Density Point

Finds an approximate mode, used to split integrals, to scale
root-finding brackets and to recenter the ratio-of-uniforms kernel.

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

A golden-section search on the compactified scale is not accurate enough
here: its tolerance is expressed in the compactified variable, whose
derivative under the tangent map is of order \\y^2\\, so a fixed
tolerance in \\t\\ corresponds to an error in \\y\\ that grows with the
location of the mode. The grid refinement instead stops on the width of
the bracket measured in \\y\\.

## See also

[`find_lp_anchor()`](https://statmodels7.github.io/distributions7/reference/find_lp_anchor.md),
[`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
