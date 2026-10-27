# Fisher Scoring, With Its Own Settings

Returns a specification of Fisher scoring for
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)'s
`method` argument. The specification records how the expected
information is obtained, and optionally a stopping rule and an iteration
limit of its own. It is
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)'s
default method.

Fisher scoring is Newton's method in which the observed Hessian is
replaced by minus the expected information, so it needs no
implementation of its own and is not an optimizer. It needs one thing
for which an optimizer has no slot, the statement of how that matrix is
obtained, and this object holds that statement.

## Usage

``` r
fisher_scoring(method = NULL, criterion = NULL, maxit = NULL)
```

## Arguments

- method:

  `NULL` (the default) or an object that selects how the expected
  information is obtained:

  - `NULL`: the closed form of the family where it has one, and
    [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
    otherwise;

  - [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md):
    the outer product of the observed scores, which costs one gradient;

  - [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md):
    quadrature (or summation over the support) of the observed Hessian
    with `identity = "hessian"`, or of minus the outer product of the
    scores with `identity = "score"`;

  - [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md):
    Monte Carlo, with the sample size `nsim` and the identity set in
    that constructor.

  Any other value, a character string included, signals an error.
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  compares the methods.

- criterion:

  A stopping rule from optimizers7, or `NULL` (the default), in which
  case
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  stops on
  [`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)
  at its default tolerance. An object that is not a `criterion` signals
  an error.

- maxit:

  An iteration limit, a single number at least 1, or `NULL` (the
  default) for the limit of
  [`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html).
  When the limit is reached the fit reports `converged = FALSE` and
  keeps the point.

## Value

An S7 object of class
[`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md),
with properties `method` (`NULL` or an object inheriting from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)),
`criterion` and `maxit`.

## Values accepted by the method argument

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
has one argument that selects how to optimize, and that argument accepts
an optimizer of optimizers7 or this object:

|  |  |
|----|----|
| `method = fisher_scoring()` | Newton's method with the expected information |
| `method = optimizers7::newton()` | Newton's method with the observed Hessian |
| `method = optimizers7::lbfgs()` | whatever that optimizer does |

## The step

Writing \\s(\eta) = \partial \ell / \partial \eta\\ for the score on the
unconstrained scale and \\\mathcal{I}(\eta) = \mathbb{E}\[-\partial^{2}
\ell / \partial \eta \partial \eta'\]\\ for the expected information
there, the step is \$\$\eta^{(t+1)} = \eta^{(t)} +
\mathcal{I}(\eta^{(t)})^{-1} s(\eta^{(t)}).\$\$ Under the diagonal
reparametrization \\\theta_i = h_i(\eta_i)\\ the information transforms
by congruence, \\\mathcal{I}(\eta) = D\\\mathcal{I}(\theta)\\D\\ with
\\D = \operatorname{diag}(h_i'(\eta_i))\\, and the first-order term of
the chain rule vanishes in expectation because the score has mean zero.
The matrix inverted is therefore positive definite wherever
\\\mathcal{I}(\theta)\\ is, and the step is a direction of ascent. A fit
therefore succeeds on a family whose observed Hessian is indefinite; the
Laplace is the standard case.

## Methods for the expected information

With `method = NULL`, the default, the fit takes the closed form of the
expected information where the family has one, and the outer product of
the observed scores
([`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md))
otherwise, as for a family defined by its density alone or a truncated
family. An object selects the method for every family, including a
family that has the closed form: the fit then uses the expectation as
that method computes it, and the closed form is not read. The methods
differ in cost and in variance:
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
costs one gradient,
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
is deterministic and dearer, and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
adds simulation noise that decreases with its `nsim`.
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
compares them.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which consumes this;
[`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md)
for the class;
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what each `method` computes;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the matrix itself;
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html),
whose iteration limit stands where this sets none.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 100, list(mu = 1, sigma = 2))

# The default and the explicit specification give the same fit.
coef(fit_distrib(d, y))
#>       mu    sigma 
#> 1.217775 1.787394 
coef(fit_distrib(d, y, method = fisher_scoring()))
#>       mu    sigma 
#> 1.217775 1.787394 

# The step rests on a congruence: the information on the link scale is
# D I(theta) D, with D the diagonal of the inverse link's derivative.
th <- list(mu = 1, sigma = 2)
It <- distrib_expected_hessian(d, y, th)
Ie <- distrib_expected_hessian(d, y, th, scale = "link")
c(theta = sum(It$sigma_sigma), link = sum(Ie$sigma_sigma))   # ratio 2^2
#> theta  link 
#>   -50  -200 

# An object in 'method' is used also by a family that has a closed form:
# here the expected information is simulated, and the fit stays close to the
# one that uses the closed form.
set.seed(2)
coef(fit_distrib(d, y, method = fisher_scoring(method = expect_mc(nsim = 2000))))
#>       mu    sigma 
#> 1.217775 1.787394 

# The stopping rule and the budget belong to the method.
fit <- fit_distrib(d, y, method = fisher_scoring(maxit = 1))
c(iterations = fit@iterations, converged = fit@converged)
#> iterations  converged 
#>          1          0 
```
