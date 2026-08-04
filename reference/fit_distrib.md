# Maximum-Likelihood Estimation

Fits a distribution to an i.i.d. sample by maximum likelihood. The
optimisation is carried out on the **link (real) scale**, where the
parameters are unconstrained, using the analytical score and information
supplied by the distribution (`scale = "link"`). Estimates are then
mapped back to the parameter scale and reported with standard errors and
confidence intervals.

## Usage

``` r
fit_distrib(
  distrib,
  y,
  start = NULL,
  method = fisher_scoring(),
  maxit = 200,
  tol = 1e-06,
  level = 0.95,
  n_start = 5
)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- start:

  Optional named list of starting values **on the parameter scale**. If
  `NULL` (default) they come from
  [`distrib_start`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
  which lets a family compute them from the data; families that do not
  say otherwise fall back to random draws, with restarts.

- method:

  How to optimise. One argument, taking one of three things:

  - [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
    the default — Newton's method with the **expected** information in
    place of the Hessian, the object carrying how that information is to
    be obtained when the family has no closed form for it;

  - an optimiser object from optimizers7, used as given and receiving
    the analytical gradient and the **observed** Hessian, so that
    `method = lbfgs(criterion = crit_grad(1e-12))` selects both the
    algorithm and the stopping rule;

  - one of the strings `"fisher"`, `"newton"` or `"bfgs"`, kept as short
    names for the three ready-made strategies. The first two fall back
    to BFGS if they fail to converge; an optimiser the caller chose is
    never silently replaced.

- maxit:

  Maximum number of iterations. Defaults to 200. An optimiser object or
  a
  [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
  carrying its own `maxit` overrides it.

- tol:

  Convergence tolerance on the score **per observation**. Defaults to
  `1e-6`. The optimiser is handed the *mean* negative log-likelihood, so
  this is what `crit_grad(tol)` tests, and a criterion supplied through
  `method` is measured on the same scale. A method object carrying its
  own stopping rule overrides this one.

  The default is not arbitrary. A line search accepts a step only when
  the objective decreases by a definite amount, and near the maximum
  that decrease is about \\\lVert U/n \rVert^2 / (2\lambda)\\ for a
  curvature \\\lambda\\. Once it falls below the rounding of the
  objective itself, about \\\varepsilon \lvert \ell/n \rvert\\, no step
  can be verified and the search stops, so the reachable floor is near
  \\\sqrt{2 \lambda \varepsilon \lvert \ell/n \rvert}\\ — of order
  `1e-8` for an objective of order one. Measured over several families,
  methods and samples that floor is usually near `1e-15` but reaches
  `1e-8`, so the default sits about two orders above it. A tighter
  tolerance asks for accuracy the arithmetic cannot certify, and whether
  a given run reaches it then depends on the platform. Nothing
  statistical is lost: a score of `1e-6` per observation places the
  estimate within a small fraction of a standard error of the maximum.

- level:

  Confidence level for the intervals. Defaults to 0.95.

- n_start:

  How many starting values to ask
  [`distrib_start`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
  for when `start` is `NULL`. Defaults to 5. A family that returns its
  own estimate returns one and ignores this.

## Value

An object of class
[`distrib_fit`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md);
see its documentation for the available components.
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`logLik()`](https://rdrr.io/r/stats/logLik.html) methods are provided.

## Details

**Why the link scale.** Optimising \\\eta \in \mathbb{R}^p\\ rather than
the constrained \\\theta\\ removes the need for box constraints: a
variance can never become negative, a probability never leaves
\\(0,1)\\. The score and information on that scale are obtained exactly
(not numerically) through the chain rule described in
[`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

**Standard errors.** The variance-covariance matrix on the link scale is
the inverse of the information at the optimum. It is mapped to the
parameter scale by the delta method,
\$\$\widehat{\mathrm{Var}}(\hat\theta) =
J\\\widehat{\mathrm{Var}}(\hat\eta)\\J, \qquad J =
\mathrm{diag}\\\left(\frac{d g^{-1}}{d\eta}\Big\|\_{\hat\eta}\right)\$\$

**Confidence intervals.** Intervals are built symmetrically on the link
scale, \\\hat\eta \pm z\_{1-\alpha/2}\\\mathrm{se}(\hat\eta)\\, and then
mapped through \\g^{-1}\\. This guarantees that the reported limits
always respect the parameter's domain (a variance interval cannot
contain negative values), which a symmetric interval on the parameter
scale would not.

## See also

[`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md),
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.html)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
d <- gaussian_distrib()
y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
fit <- fit_distrib(d, y)
fit
coef(fit)
vcov(fit)

# a bounded parameter: the interval never leaves (0, 1)
b <- bernoulli_distrib()
fit_distrib(b, rbinom(50, 1, 0.9))
} # }
```
