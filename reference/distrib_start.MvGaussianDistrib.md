# The Maximum Likelihood Estimate as a Starting Value

Returns the sample mean and the sample covariance, carried onto the
distribution's matrix parametrization. For an unstructured covariance
those **are** the maximum likelihood estimate, so the fit begins at the
answer and confirms it in one step; for a structured one they are the
closest matrix the parametrization can represent, which is a great deal
nearer than the origin.

One starting value is returned whatever `n_start` asks for. There is
nothing to explore when the first point is the estimate.

## Arguments

- distrib:

  An
  [`MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  The response, an \\n \times p\\ matrix.

- n_start:

  Ignored: one starting value is returned, and it is the estimate.

- ...:

  Unused.

## Value

A list of length 1 holding one named parameter list: the \\p\\ location
components followed by the structure's free values, named and ordered as
`distrib@params`.

## Details

Where the object parametrizes the **precision**, the sample covariance
is inverted first: the matrix the structure has to represent is
\\\hat\Sigma^{-1}\\, not \\\hat\Sigma\\.

The covariance divides by \\n\\, which is the maximum likelihood
estimator, and its eigenvalues are floored by
[`mv_moment_start()`](https://statmodels7.github.io/distributions7/reference/mv_moment_start.md)
so that a singular sample covariance still produces a usable point.
Carrying it onto the structure goes through
[`param_free_or_fit()`](https://statmodels7.github.io/distributions7/reference/param_free_or_fit.md),
which is exact where
[`parameters7::param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.html)
succeeds and a least-squares fit where it does not.

## See also

[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic;
[`distrib_start.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_start.MvStudentTDistrib.md),
which starts from this;
[`mv_moment_start()`](https://statmodels7.github.io/distributions7/reference/mv_moment_start.md)
and
[`param_free_or_fit()`](https://statmodels7.github.io/distributions7/reference/param_free_or_fit.md)
for the two steps.
