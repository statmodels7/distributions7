# Default Numerical Quantile Function for Continuous Distributions

The fallback for a continuous family that implements no analytical
quantile function: it inverts
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
by root finding, which on a family with no analytical distribution
function means inverting the quadrature of
[`distrib_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md).

## Arguments

- distrib:

  An object inheriting from `continuous_distrib` that registers no
  method of its own.

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  under `log.p`. The bounds of the support are returned at 0 and 1, and
  `NaN` for a value outside \\\[0, 1\]\\.

- theta:

  A named list of parameters, each of length 1 or `length(p)`.

- lower.tail:

  Logical of length 1. `TRUE`, the default, treats `p` as \\P(Y \le
  q)\\; `FALSE` treats it as the upper-tail probability.

- log.p:

  Logical of length 1. `TRUE` treats `p` as a logarithm.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, the length of the recycled `p` and
`theta`.

## Details

Every probability is inverted at once. The bracket starts at an
approximate mode, found by the grid search of
[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md)
run on every distinct parameter setting in one evaluation of the density
per refinement
([`quantile_anchors()`](https://statmodels7.github.io/distributions7/reference/quantile_anchors.md)),
and expands geometrically, with a first step equal to the reciprocal of
the density there, so that a sharply peaked family and a diffuse one
take a comparable number of expansions; after thirty doublings the step
is multiplied by a thousand, so that a quantile far in a tail is
reached. The expansion stops at a quarter of the largest double: a root
beyond it is returned as `Inf` (or `-Inf`), the quadrature of a
distribution function being unable to reach the mass past the largest
double.

The probability is first turned into the logarithms of both tails by
[`quantile_log_tails()`](https://statmodels7.github.io/distributions7/reference/tail_helpers.md),
without forming \\1 - p\\ or exponentiating a logarithm, and each root
is solved on the tail with the smaller probability: \\\log F(q) = \log
P\\ below the median and \\\log S(q) = \log P\\ above it, with \\S\\ the
upper tail that the family's
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
returns. Inside the bracket a Newton step on that equation, whose
derivative is \\\pm f/P\\, is taken where it stays inside the bracket
and a bisection where it does not, so every iteration shrinks the
bracket or converges; the tolerance is relative on the log scale, so a
probability far below machine epsilon is inverted to full relative
accuracy. Each iteration is one call of
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
and one of
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
over every probability not yet inverted.

The root finding is applied to a distribution function that is itself
computed by quadrature, so the quantile is less accurate than the
distribution function that it inverts.

## See also

[`distrib_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md),
the function it inverts;
[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md)
for the bracket's starting point;
[`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md),
which tests whether this method is what a caller would reach;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.
