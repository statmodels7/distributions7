# Maximum-Likelihood Estimation for a Distribution

Fits a `distrib` object to an i.i.d. sample by maximum likelihood and
returns a
[`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
carrying the estimates, their standard errors and confidence limits on
both scales. The optimization runs on the link scale, where the
parameters are unconstrained, and it is driven by the analytical score
and information that the family supplies at `scale = "link"`; the
estimates are then mapped back and reported on the parameter scale.

On a Gaussian the answer is the closed-form estimate to the printed
digit: \\\hat\mu = \bar y\\, \\\hat\sigma^2 = \frac{1}{n}\sum (y_i -
\bar y)^2\\, with \\\mathrm{se}(\hat\mu) = \hat\sigma/\sqrt{n}\\ and
\\\mathrm{se}(\hat\sigma) = \hat\sigma/\sqrt{2n}\\.

## Usage

``` r
fit_distrib(
  distrib,
  y,
  start = NULL,
  method = fisher_scoring(),
  level = 0.95,
  n_start = 5,
  threads = numericals7::n_threads()
)
```

## Arguments

- distrib:

  An object inheriting from `distrib`: a univariate family, a
  multivariate one, or any wrapper of either.

- y:

  A numeric vector of observations, or an \\n \times p\\ matrix for a
  multivariate family. Every observation is read at the same \\\theta\\.

- start:

  Optional named list of starting values on the parameter scale. `NULL`,
  the default, takes the starting values from
  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
  which computes them from the data for a family that registers a method
  and otherwise draws them from the parameter domains, with the restarts
  below. A value that is neither `NULL` nor a list signals an error
  naming the argument, because `start` sits before `method` in the
  signature and an optimizer passed positionally lands here.

- method:

  An object that selects how to optimize. The argument takes one of two
  kinds of object:

  - [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
    the default: Newton's method with the expected information in place
    of the observed Hessian, the object carrying the method
    ([ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md))
    that obtains that information. With `method = NULL` the closed form
    of the family is used where there is one, and
    [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
    otherwise; a method given there is used for every family. Fisher
    scoring falls back to BFGS if it fails to converge;

  - an optimizer object from optimizers7, used as given and receiving
    the analytical gradient and the observed Hessian, so that
    `method = optimizers7::lbfgs(criterion = optimizers7::crit_grad(1e-12))`
    selects the algorithm and the stopping rule together. A stopping
    rule that the optimizer cannot evaluate is rejected here, where the
    message can name it. An optimizer that the caller chose is never
    replaced.

  Any other value, a character string included, signals an error.

  The iteration limit and the stopping rule belong to the method and are
  set there, on an optimizer through its own `maxit` and `criterion` and
  on
  [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
  through the same two arguments. For
  [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
  the rule is
  [`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)
  at its default tolerance when the method does not set a rule; an
  optimizer object keeps its own default.

- level:

  Confidence level for `ci` and `ci_eta`, a single number in \\(0, 1)\\;
  any other value signals an error. Defaults to 0.95. Any other level is
  available afterwards from
  [`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
  without refitting.

- n_start:

  How many starting values to request from
  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
  when `start` is `NULL`. Defaults to 5. A family that returns its own
  estimate returns one and ignores this. It is ignored entirely when
  `start` is given.

- threads:

  How many threads the fit may use, as
  [`numericals7::n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.html)
  constructs it. The default, `n_threads(1)`, is sequential and takes
  the sequential code path. The count reaches the family's compiled
  per-observation kernels as an argument; the result does not depend on
  it, bit for bit, because every parallel region decomposes its work
  over the elements of its output and never splits a reduction.

## Value

An S7 object of class
[`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md).
Its `coefficients`, `se` and `ci` are on the parameter scale and its
`eta`, `se_eta` and `ci_eta` on the link scale; `loglik`, `aic` and
`bic` are computed from the unscaled log-likelihood at the estimate, and
`converged`, `score`, `iterations`, `method`, `criterion` and `counts`
record what the optimizer did. See that page for every component.

Signals an error when every starting value fails to produce a usable
run, and a separate one when every run of a discrete family reached a
positive log-likelihood, which is impossible for a product of
probabilities and indicates that the mass function has broken down at
the parameters reached.

## Why the link scale

Optimizing \\\eta \in \mathbb{R}^p\\ in place of the constrained
\\\theta\\ removes the need for box constraints: a scale cannot become
negative and a probability cannot leave \\(0,1)\\, because every link
maps onto the interior of its parameter's domain. The score and the
information on that scale follow from those on the parameter scale
through the chain rule of
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md),
and they are exact wherever the family's own derivatives are.

## Scaling of the objective

The optimizer receives \\-\ell(\eta)/n\\, and the gradient and the
Hessian are divided by \\n\\ as well. Scaling by a positive constant
moves neither the maximum nor any Newton step, since \\H^{-1}g\\ is
unchanged when both are divided by \\n\\. The scaling changes the
meaning of a threshold: a tolerance on the gradient of this objective is
a tolerance on the score per observation whatever the sample size, so
the same rule means the same thing at \\n = 10\\ and at \\n = 10^7\\.
`loglik`, the information and every standard error are recomputed
unscaled at the optimum.

## Standard errors and intervals

The variance matrix on the link scale is the inverse of the observed
information at \\\hat\eta\\, the negative Hessian of the log-likelihood,
for every family and every `method`. The matrix with which the optimizer
ended is not used, and another matrix is selected after the fit with the
`information` argument of
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
and
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
(see
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)).
The delta method carries the matrix to the parameter scale,
\$\$\widehat{\mathrm{Var}}(\hat\theta) =
J\\\widehat{\mathrm{Var}}(\hat\eta)\\J, \qquad J =
\mathrm{diag}\\\left(\frac{dg^{-1}}{d\eta}\Big\|\_{\hat\eta}\right).\$\$
Intervals are built symmetrically on the link scale, \\\hat\eta \pm
z\_{1-\alpha/2}\\\mathrm{se}(\hat\eta)\\, and mapped through \\g^{-1}\\,
sorting the pair in case the link decreases. The limits therefore
respect the parameter's domain, which a symmetric interval on the
parameter scale would not.

## Restarts, the fallback and the tie-break

Each starting value is tried in turn, and the search stops once two runs
that converged have reached the best objective found so far, to within
\\10^{-8}\\ times the larger of 1 and the absolute value of the mean
negative log-likelihood. A single converged run is not taken as the
answer because it may be a degenerate stationary point. Where fewer than
two runs converge every start is tried and the best is kept. Fisher
scoring falls back to BFGS from the same starting value when it fails;
an optimizer object that the caller supplied does not. Among the runs
that finish, a run that converged beats one that did not converge and
the objective breaks ties, so the fit reports the best run and not the
last. A run of a discrete family whose log-likelihood came back positive
is discarded before any comparison.

## See also

[`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
for the object returned;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for the default method;
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for where the starting values come from;
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
for the chain rule that the score uses;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family before fitting it;
[`optimizers7::minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.html),
which runs the search.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
fit <- fit_distrib(d, y)
fit
#> Maximum-likelihood fit: gaussian1
#> Observations: 500   Log-likelihood: -1264   AIC: 2532   BIC: 2541
#> Method: Fisher scoring   iterations: 2   evaluations: f 3, g 3   time: 6 ms
#> Converged: yes (gradient (max-norm) < 1e-06)
#> 
#> Parameter scale:
#>       Estimate Std. Error   2.5%  97.5%
#> mu      2.0679     0.1356 1.8021 2.3338
#> sigma   3.0327     0.0959 2.8505 3.2267
#> 
#> Link scale:
#>       Estimate Std. Error   2.5%  97.5%
#> mu      2.0679     0.1356 1.8021 2.3338
#> sigma   1.1095     0.0316 1.0475 1.1714

# The Gaussian MLE is closed form, and the fit reaches it.
all.equal(coef(fit)[["mu"]], mean(y))
#> [1] TRUE
all.equal(coef(fit)[["sigma"]], sqrt(mean((y - mean(y))^2)))
#> [1] TRUE

# So are the standard errors: sigma/sqrt(n) and sigma/sqrt(2n).
s <- coef(fit)[["sigma"]]
all.equal(unname(fit@se), c(s / sqrt(500), s / sqrt(2 * 500)))
#> [1] TRUE

# A bounded parameter: the interval is built on the link scale and mapped
# back, so it cannot contain a probability outside (0, 1).
b <- bernoulli_distrib()
fb <- fit_distrib(b, rbinom(50, 1, 0.9))
rbind(link = confint(fb, scale = "link"), parameter = confint(fb))
#>         2.5%     97.5%
#> mu 1.2732887 3.1211605
#> mu 0.7813052 0.9577572

# A non-regular family. The Laplace's observed curvature in the location is
# zero almost everywhere, so Newton's method has nothing to invert; Fisher
# scoring uses the information instead and reaches the closed-form estimates.
yl <- distrib_rng(laplace_distrib(), 400, list(mu = 0, sigma = 1))
fl <- fit_distrib(laplace_distrib(), yl)
c(fitted = coef(fl)[["mu"]], median = median(yl))
#>    fitted    median 
#> -0.120917 -0.120917 
c(fitted = coef(fl)[["sigma"]], mad = mean(abs(yl - median(yl))))
#>   fitted      mad 
#> 1.001251 1.001251 
```
