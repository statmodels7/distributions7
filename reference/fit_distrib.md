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
  method = c("fisher", "newton", "bfgs"),
  maxit = 200,
  tol = 1e-10,
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
  `NULL` (default) starting values are drawn with
  [`generate_random_theta`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md),
  with a few random restarts on failure; supplying a sensible `start`
  makes convergence faster and more reliable.

- method:

  Optimisation method. `"fisher"` (default) uses Fisher scoring with the
  expected information, `"newton"` uses the observed Hessian, and
  `"bfgs"` uses [`optim`](https://rdrr.io/r/stats/optim.html) with the
  analytical gradient. Fisher scoring and Newton fall back to BFGS if
  they fail to converge.

- maxit:

  Maximum number of iterations. Defaults to 200.

- tol:

  Convergence tolerance on the score and on the log-likelihood
  increment. Defaults to `1e-10`.

- level:

  Confidence level for the intervals. Defaults to 0.95.

- n_start:

  Number of random restarts attempted when `start` is `NULL` and the
  first attempt fails. Defaults to 5.

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
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)

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
