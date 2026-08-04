# The Step a Skew t Differences the Degrees of Freedom With

Returns the finite-difference step used for the derivatives in \\\nu\\,
relative to \\\nu\\ itself and floored so that it stays meaningful for a
small number of degrees of freedom.

## Usage

``` r
skewt_nu_step(nu)
```

## Arguments

- nu:

  The degrees of freedom.

## Value

A numeric vector.

## Details

The relative step \\10^{-3}\\ is measured rather than assumed. Swept
over \\\nu\\ from 2 to 30 and sample sizes from 500 to 4000, it is where
the truncation error of the five-point stencil of
[`fd5_first`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
has fallen to the level of the rounding error and the two are balanced;
a smaller step is dominated by rounding, which the stencil amplifies by
\\18/(12h)\\, and a larger one by truncation, which grows as \\h^4\\.
