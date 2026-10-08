# The Pseudo-Huber Starts at a Finite Shape Below Gaussian Kurtosis

Returns \\\mu = \bar y\\, \\\nu = 10\\ and the \\\sigma\\ that gives the
sample variance at that \\\nu\\, when the sample excess kurtosis is
below 0.05, and an empty list otherwise.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md),
  or a `PseudoHuber2Distrib` object, from
  [`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md).

- y:

  The response, a numeric vector.

- ...:

  Unused.

## Value

A named list with `mu`, `sigma` and `nu`, or an empty list.

## Details

The excess kurtosis of the pseudo-Huber family lies in \\(0, 3)\\, and
it tends to zero in the gaussian limit \\\nu \to \infty\\, where the
variance \\\sigma^2 \sqrt{\nu} K_2(\sqrt{\nu}) / K_1(\sqrt{\nu})\\ stays
finite only if \\\sigma \to 0\\. When the sample excess kurtosis of the
response is zero or negative, the intercept-only maximum likelihood fit
reaches that limit, and a regression started there stays on the ridge
\\\sigma \propto \nu^{-1/4}\\. Measured on ten simulated regressions
with a strong covariate: the intercept-only fit ended at \\\log\nu\\
between 27 and 40 in the six samples with a non-positive excess kurtosis
and at most 6.0 in the others, and the regression started there stopped
7.3 to 15.8 log-likelihood units below its maximum in four of the six.
From \\\nu = 1\\, 3 or 10 every one of the ten reached its maximum.
\\\nu = 10\\ gives an excess kurtosis of 0.85.

The same method is registered on
[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md),
where the returned `sigma` is the sample standard deviation.

## See also

[`distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.md)
for the generic,
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).
