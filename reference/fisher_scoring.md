# Fisher Scoring, With Its Own Settings

Returns a specification of Fisher scoring for
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)'s
`method` argument, carrying how the expected information is to be
obtained when the family has no closed form for it, and optionally a
stopping rule and an iteration limit of its own. It is
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)'s
default method.

Fisher scoring is Newton's method with the observed Hessian replaced by
minus the expected information, so it needs no implementation of its own
and is not an optimizer. What it does need, and an optimizer has no slot
for, is the statement of how that matrix is to be obtained, and that is
what this object holds.

## Usage

``` r
fisher_scoring(
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000,
  criterion = NULL,
  maxit = NULL
)
```

## Arguments

- approx:

  How the expectation is approximated when the family has no closed-form
  expected information. `"bartlett"` (the default) is the outer product
  of the score, for which `"opg"` is an accepted synonym returning the
  identical value; `"integrate"` is quadrature of the observed
  information; `"mc"` is Monte Carlo. Matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html), so
  anything else signals an error listing the four. See
  [`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- nsim:

  Monte Carlo sample size, a single positive finite number. Read only
  when `approx = "mc"`. Defaults to 10000. A vector, a non-finite value
  or a number below 1 signals an error.

- criterion:

  A stopping rule from optimizers7, or `NULL` (the default) to leave
  [`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)'s
  own in force. Anything that is not a `criterion` object signals an
  error.

- maxit:

  An iteration limit, a single number at least 1, or `NULL` (the
  default) for the same. When the limit is reached the fit reports
  `converged = FALSE` and keeps the point.

## Value

An S7 object of class
[`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md),
with properties `approx`, `nsim`, `criterion` and `maxit`.

## One argument, three kinds of thing

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
takes one argument saying how to optimize, and it takes an optimizer of
optimizers7 or this:

|  |  |
|----|----|
| `method = fisher_scoring()` | Newton's method with the **expected** information |
| `method = optimizers7::newton()` | Newton's method with the **observed** Hessian |
| `method = optimizers7::lbfgs()` | whatever that optimizer does |

## The step

Writing \\s(\eta) = \partial \ell / \partial \eta\\ for the score on the
unconstrained scale and \\\mathcal{I}(\eta) = \mathbb{E}\[-\partial^{2}
\ell / \partial \eta \partial \eta'\]\\ for the expected information
there, \$\$\eta^{(t+1)} = \eta^{(t)} + \mathcal{I}(\eta^{(t)})^{-1}
s(\eta^{(t)}).\$\$ Under the diagonal reparametrization \\\theta_i =
h_i(\eta_i)\\ the information transforms by congruence,
\\\mathcal{I}(\eta) = D\\\mathcal{I}(\theta)\\D\\ with \\D =
\operatorname{diag}(h_i'(\eta_i))\\, the first-order term of the chain
rule vanishing in expectation because the score has mean zero. The
matrix inverted is therefore positive definite wherever
\\\mathcal{I}(\theta)\\ is, and the step is a direction of ascent. A fit
therefore succeeds on a family whose observed Hessian is indefinite, the
Laplace being the standard case.

## Which strategy, and what it costs

`approx` is read only where the family has no closed-form expected
information; nine families in this package are in that position, and for
them the choice is a trade of cost against noise. Measured on a skew
normal at \\(\mu, \sigma, \alpha) = (0, 1, 3)\\ over 200 observations,
`"bartlett"` and `"integrate"` agree to \\8\times 10^{-14}\\ on every
component at 0.28 s and 0.18 s, and `"mc"` at `nsim = 2000` lands within
about one per cent of them for no measurable time.

A family that supplies the expectation in closed form ignores `approx`
entirely, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
signals an error when one is given there, so that a caller cannot come
to believe a fit used a strategy it never read.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which consumes this;
[`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md)
for the class;
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what each `approx` computes;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the matrix itself;
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html),
whose defaults stand where this sets none.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 100, list(mu = 1, sigma = 2))

# The default, and the same thing said explicitly.
coef(fit_distrib(d, y))
#>       mu    sigma 
#> 1.217775 1.787394 
coef(fit_distrib(d, y, method = fisher_scoring()))
#>       mu    sigma 
#> 1.217775 1.787394 

# The congruence the step rests on: the information on the link scale is
# D I(theta) D, with D the diagonal of the inverse link's derivative.
th <- list(mu = 1, sigma = 2)
It <- distrib_expected_hessian(d, y, th)
Ie <- distrib_expected_hessian(d, y, th, scale = "link")
c(theta = sum(It$sigma_sigma), link = sum(Ie$sigma_sigma))   # ratio 2^2
#> theta  link 
#>   -50  -200 

# 'approx' is read only where the family has no closed form. On a skew
# normal the two deterministic strategies agree; Monte Carlo is noisy.
sn <- skewnormal1_distrib()
set.seed(2)
ys <- distrib_rng(sn, 40, list(mu = 0, sigma = 1, alpha = 3))
ths <- list(mu = 0, sigma = 1, alpha = 3)
vapply(c("bartlett", "integrate"), function(a) {
  sum(distrib_expected_hessian(sn, ys, ths, approx = a)$alpha_alpha)
}, numeric(1))
#>  bartlett integrate 
#> -1.085353 -1.085353 

# The same argument on a family that HAS one is rejected, not ignored.
try(fit_distrib(d, y, method = fisher_scoring(approx = "mc")))
#> Error : 'gaussian1' computes its expected information in closed form, so the 'approx'
#>   of fisher_scoring() would be ignored. Use fisher_scoring() with no
#>   arguments: the fit will take the exact expression.

# The stopping rule and the budget belong to the method.
fit <- fit_distrib(d, y, method = fisher_scoring(maxit = 1))
c(iterations = fit@iterations, converged = fit@converged)
#> iterations  converged 
#>          1          0 
```
