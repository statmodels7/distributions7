# Score of a Truncated Distribution

The parent's score recentred by its truncated mean, \\d_i \ell_T =
s_i(y) - m_i\\.

## Usage

``` r
trunc_gradient(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A truncated distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list, one component per parameter.

## Details

Truncation adds no parameter – the endpoints are constants, like a
binomial's size – but it does add a \\\theta\\-dependent normalising
constant, and the recentring is that constant's contribution. The
support does not depend on \\\theta\\, which is what licenses
differentiating \\Z\\ under the integral sign and keeps truncation at
fixed points a regular problem.
