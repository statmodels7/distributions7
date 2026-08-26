# Mean of a von Mises in the Resultant Length

Returns the ordinary expectation of \\Y\\ as a number on \\\[-\pi,
\pi)\\, obtained numerically by delegating to
[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
at the implied concentration.

It is **not** a circular quantity, and neither parameter describes it.
\\\mu\\ is the mean *direction* and \\\rho\\ the mean resultant length;
\\\mathbb{E}\[Y\]\\ differs from \\\mu\\ whenever \\\mu \ne 0\\, because
the interval is cut at \\\pm\pi\\ and the density is not symmetric about
\\\mu\\ on it. Compute the circular mean of a sample as
`atan2(mean(sin(z)), mean(cos(z)))`, which recovers \\\mu\\.

## Arguments

- x:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).
  The argument is named `x` because the generic is
  [`base::mean()`](https://rdrr.io/r/base/mean.html).

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1. `rho` must lie in \\(0, 1)\\. Aligned and validated by name,
  so a missing or out-of-bounds component throws.

- ...:

  Passed on to
  [`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
  which reads the quadrature's settings.

## Value

A numeric vector of length 1, the ordinary mean of \\Y\\.

## See also

[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
which this calls;
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md),
whose page discusses the same distinction; and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the other ordinary moments.
