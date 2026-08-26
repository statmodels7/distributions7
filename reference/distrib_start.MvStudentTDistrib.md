# The Gaussian Estimate as a Starting Value for a t

Returns the sample mean and the sample covariance with the degrees of
freedom set at 8, which is heavy-tailed enough to be worth fitting and
light enough for the second moment to exist. The Gaussian estimate is
this family's limit as \\\nu\\ grows, so it is where a finite \\\nu\\ is
looked for from.

One starting value is returned whatever `n_start` asks for.

## Arguments

- distrib:

  An
  [`MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  The response, an \\n \times p\\ matrix.

- n_start:

  Ignored: one starting value is returned.

- ...:

  Unused.

## Value

A list of length 1 holding one named parameter list: the \\p\\ location
components, the structure's free values, then `nu` at 8, named and
ordered as `distrib@params`.

## Details

The scale matrix is not the covariance. For a \\t\\ with \\\nu\\ degrees
of freedom \\\mathrm{Var}(Y) = \nu\Sigma/(\nu-2)\\, so the sample
covariance is multiplied by \\(\nu_0-2)/\nu_0 = 0.75\\ before it is
carried onto the structure. A starting value that confused the two would
begin with a scale a third too large.

## See also

[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic;
[`distrib_start.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_start.MvGaussianDistrib.md),
the limit this starts from;
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the scale matrix and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the covariance.
