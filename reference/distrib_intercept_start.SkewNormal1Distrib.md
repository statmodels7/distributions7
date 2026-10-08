# The Direct Skew Normal Starts at Its Moment Estimate

Returns the moment estimate of \\(\mu, \sigma, \alpha)\\, with the
sample skewness held at 0.9, when the absolute sample skewness is at
least 0.9, and an empty list otherwise.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  The response, a numeric vector.

- ...:

  Unused.

## Value

A named list with `mu`, `sigma` and `alpha`, or an empty list.

## Details

The skewness of the skew normal is bounded by \\\gamma\_{\max} \approx
0.9953\\. When the sample skewness of the response is beyond that bound,
the intercept-only maximum likelihood fit has \\\alpha \to \infty\\ and
\\\mu\\ at the smallest observation, which is the half-normal limit. In
the direct parametrization the curvature in \\\mu\\ there grows as
\\\alpha^2\\, so a regression started from that point takes steps of
order \\1/\alpha^2\\ and does not reach its maximum. The moment estimate
with the skewness held at 0.9 gives \\\alpha = 6.3\\ and the matching
location and scale. Measured on
[`MASS::Cars93`](https://rdrr.io/pkg/MASS/man/Cars93.html), price on
horse power (marginal skewness 1.48): the intercept-only start was at
\\\alpha = 1.4 \times 10^6\\ and the REML criterion was unavailable at
the starting point, and from the moment estimate the fit converged to
the marginal log-likelihood of
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
on the same model.

## See also

[`distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.md)
for the generic,
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).
