# Default Numerical Quantile Function for Continuous Distributions

The fallback for a continuous family that implements no analytical
quantile function: it inverts
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
by root-finding, which on a family with no analytical distribution
function either means inverting the quadrature of
[`distrib_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md).
Measured on a Gamma defined by its density alone, it agrees with
[`stats::qgamma()`](https://rdrr.io/r/stats/GammaDist.html) to between
\\4.1\times10^{-11}\\ and \\4.5\times10^{-10}\\ relative, and the round
trip \\F(F^{-1}(p)) - p\\ closes to \\10^{-11}\\.

## Arguments

- distrib:

  An object inheriting from `continuous_distrib` that registers no
  method of its own.

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  under `log.p`. The bounds of the support are returned at 0 and 1.

- theta:

  A named list of parameters, each of length 1 or `length(p)`.

- lower.tail:

  Logical of length 1. `TRUE`, the default, treats `p` as \\P(Y \le
  q)\\; `FALSE` treats it as the upper tail and inverts \\1 - p\\.

- log.p:

  Logical of length 1. `TRUE` treats `p` as a logarithm and
  exponentiates it first, so a probability that has underflowed is not
  recovered.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, the length of the recycled `p` and
`theta`.

## Details

The bracket starts at an approximate mode from
[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md)
and expands geometrically, with the step scaled by the density height
there so that a sharply peaked family and a diffuse one take a
comparable number of expansions. The mode and its scale are computed
once per distinct parameter setting and reused across the probabilities
that share it.

Two layers of numerical work stack here, which is why the accuracy is
four orders coarser than the distribution function's: the root-finder
can only be as accurate as the function it inverts.

## See also

[`distrib_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md),
the function it inverts;
[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md)
for the bracket's starting point;
[`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md),
which asks whether this method is what a caller would reach;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.
