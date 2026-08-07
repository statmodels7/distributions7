# Locate a High-Density Point of a Bare Log-Density

The search behind
[`find_pdf_anchor`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md),
expressed on a plain log-density over an interval.

## Usage

``` r
find_lp_anchor(lp_raw, b)
```

## Arguments

- lp_raw:

  A function giving the log-density at a numeric vector.

- b:

  A length-2 numeric vector, the interval to search.

## Value

A single number.

## Details

Kept separate from the `distrib` object so that it can also be applied
to a reparametrized density, which is what the divergence-removing
transforms in this file produce and which has no distribution object of
its own.

## See also

[`find_pdf_anchor`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md)
