# Default Numerical Hessian for `distrib` Objects

The fallback for a family that implements no analytical Hessian: second
differences of `distrib_pdf(..., log = TRUE)` through
[`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md).
A diagonal component takes the three-point stencil
\\\[\ell(\theta_i+h) - 2\ell(\theta_i) + \ell(\theta_i-h)\]/h^2\\ and an
off-diagonal one the four-point mixed stencil, so both are a **single**
difference of the log-density and neither is a difference of the
gradient.

## Arguments

- distrib:

  An object inheriting from `distrib` that registers no method of its
  own.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  Handled by the generic after dispatch; this method always returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list of Hessian component vectors, each of length `length(y)`,
keyed by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
which puts the diagonal first and is **not** the lexicographic keying
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
uses above order 2.

## The step and the cost

The step is \\h = \varepsilon^{1/4}\max(1, \|\theta_i\|) \approx
1.22\times10^{-4}\\, twenty times the gradient's: a second difference
divides by \\h^2\\, so rounding grows as \\1/h^2\\ and the optimum moves
out.
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
applies the same boundary clamp. One Hessian costs \\2p\\ evaluations
for the diagonal and \\4\\ per distinct pair, which is 6 in all for a
two-parameter family.

## What it delivers

Measured on a Gamma in its mean and dispersion at \\(\mu, \sigma^2) =
(2, 0.7)\\, against the family's own closed form: the three components
agree to \\2.3\times10^{-9}\\, \\2.7\times10^{-8}\\ and
\\3.6\times10^{-8}\\ relative, two to three digits worse than the
gradient's. That is the price of a second difference.

## See also

[`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md),
which does the differencing;
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
for the boundary rule;
[`distrib_gradient.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.distrib.md)
for the order below;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the expectation of this.
