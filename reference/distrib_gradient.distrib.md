# Default Numerical Gradient for `distrib` Objects

The fallback for a family that implements no analytical score: the
gradient of `distrib_pdf(..., log = TRUE)` by one **central difference**
per parameter, through
[`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md).
This is why
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
is the only compulsory method of the package: a family that defines the
density alone gets a score, an information, four orders of derivative
and a fit.

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

A named list with one numeric vector per parameter, keyed by
`distrib@params`, each of length `length(y)`.

## The stencil, the step and the cost

Each component is \\\[\ell(\theta_i + h) - \ell(\theta_i - h)\]/(2h)\\,
so one gradient costs \\2p\\ evaluations of the log-density. The step is
\\h = \varepsilon^{1/3}\max(1, \|\theta_i\|) \approx 6.06\times10^{-6}\\
at a parameter of order one, which balances the \\O(h^2)\\ truncation of
a central difference against a rounding term growing as \\1/h\\.
Parameter domains here are open and a step through zero returns `NaN`
from the density for reasons that look like a defect in the family, so
near a finite boundary the step is kept inside:
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
offers it cut to 49\\ or scaled on that distance, and
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md)
keeps whichever agrees better with itself at half its own step. On a
gamma at a dispersion of \\10^{-6}\\ that choice is worth five orders,
the relative error going from 3.5e-01 to 1.1e-06.

## What it delivers

Measured on a Gamma in its mean and dispersion at \\(\mu, \sigma^2) =
(2, 0.7)\\, against the family's own closed form: the two components
agree to \\1.3\times10^{-11}\\ and \\9.7\times10^{-11}\\ relative, which
is the \\O(h^2)\\ the step promises.

## See also

[`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md),
which does the differencing;
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
for the boundary rule;
[`distrib_hessian.distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.distrib.md)
for the order above;
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.
