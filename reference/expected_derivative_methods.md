# Methods for the Expectation of a Derivative

An expected derivative of the log-density, such as the expected Hessian
(the negative of the Fisher information) or the expected third
derivative, is the expectation \\\mathbb{E}\[\\\cdot\\\]\\ of a function
of the response taken under the distribution at the parameters. A family
may write that expectation in closed form. When it does not, or when a
different route is wanted, the expectation is obtained by one of three
methods, which are objects:
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md).
All of them inherit from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
and are passed through the argument `method` of
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
the derivatives of the expected information
([`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and the functions of its family),
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md).

The argument `method` has one rule. If it is `NULL`, which is the
default, the function returns the closed form of the family where the
family has one, and otherwise the default route of that function (the
table below). If it is an object, the object is always used, also for a
family that has the expectation in closed form; the closed form is then
bypassed, and the result is the expectation as the method computes it.

## Value

Nothing. This page documents the argument `method` shared by the
functions named above; the value returned is theirs.

## Two identities and three ways

Let \\\ell\\ be the log-density and \\\ell\_{i}, \ell\_{ij}, \dots\\ its
derivatives with respect to the parameters. The expected derivative of
order \\k\\ has two expressions inside the expectation, which are equal
in expectation for a regular model.

- The Hessian identity:

  The expectation of the observed derivative of order \\k\\ itself,
  \$\$\mathbb{E}\[\ell\_{i_1 \dots i_k}\].\$\$ At order 2 this is the
  expected Hessian read literally.

- The score identity:

  The Bartlett partition identity, which writes the same expectation
  through products of derivatives of lower order. The identity states
  that \$\$\sum\_{\pi} \mathbb{E}\left\[\prod\_{B \in \pi}
  \ell_B\right\] = 0,\$\$ where \\\pi\\ runs over the partitions of the
  index set \\\\i_1,\dots,i_k\\\\ and \\\ell_B\\ is the derivative with
  respect to the parameters indexed by the block \\B\\. The single-block
  partition is the target, so that \$\$\mathbb{E}\[\ell\_{i_1 \dots
  i_k}\] = -\sum\_{\pi,\\ \|\pi\| \ge 2} \mathbb{E}\left\[\prod\_{B \in
  \pi} \ell_B\right\].\$\$ At order 2 this is \\-\mathbb{E}\[\ell_i
  \ell_j\]\\, at order 3 \\-\mathbb{E}\[\ell\_{ij}\ell_k +
  \ell\_{ik}\ell_j + \ell\_{jk}\ell_i + \ell_i\ell_j\ell_k\]\\, and at
  order 4 the sum over the partitions of four indices into at least two
  blocks. The target order is never needed.

The expectation of either quantity is then obtained in one of three
ways: computed by quadrature (or by a series over the support of a
discrete family), simulated, or estimated from the sample in hand. The
methods fill the table as follows.

|  |  |  |  |
|----|----|----|----|
|  | computed by quadrature | simulated | estimated from the sample |
| Hessian identity | `expect_integrate(identity = "hessian")` | `expect_mc(identity = "hessian")` | not available |
| score identity | `expect_integrate(identity = "score")` | `expect_mc(identity = "score")` | [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md) |

The estimate from the sample has no Hessian-identity entry because the
observed derivative at the data, averaged over the observations, is the
observed information and not an expectation. Reading the score identity
at the data gives the empirical outer product of the scores, which is
the BHHH estimator at order 2.

- [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md):

  Computes the expectation by batched adaptive quadrature for a
  continuous family and by series summation for a discrete one. It is
  deterministic, and it is normally the most accurate method. With
  `identity = "hessian"` it costs one integral per component and needs
  the derivative of the target order, which every family of the package
  supplies or approximates by finite differences. With
  `identity = "score"` it costs one integral per partition and component
  and needs only the lower orders, which is the reason to choose it
  where only those are available in closed form. The tolerances and
  limits are arguments of the constructor. For a family with several
  parameters varying by observation, one integral is computed per
  distinct parameter configuration.

- [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md):

  Draws `nsim` observations from the distribution, with
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
  at each parameter configuration and averages the quantity of the
  identity over them. It estimates the same quantity as
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  with the same `identity`. The error decreases as
  \\1/\sqrt{n\_{sim}}\\, and the result is reproducible only if the seed
  is fixed. The cost is that of
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
  times `nsim`, per parameter configuration. It is the simplest method
  where quadrature performs poorly (heavy tails, awkward supports); for
  a distribution whose generator inverts a slow quantile function,
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  is cheaper.

- [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md):

  Reads the score identity at each observation and takes no expectation:
  nothing is integrated and nothing is simulated. At order 2 the
  quantity is \\-\ell_i \ell_j\\, evaluated with one call to
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md);
  at orders 3 and 4 it is the sum of products of the identity, evaluated
  with the derivatives of orders 1 to \\k-1\\. Each observation
  contributes a noisy rank-one-like term, and the sum over \\n\\
  observations converges to \\n\\ times the expectation at the rate
  \\n^{-1/2}\\. At order 2 the matrix is positive semidefinite by
  construction, which makes it adequate for the step of Fisher scoring,
  and it has the same fixed point as the expected information. The
  variance grows with the order, because the terms are products of up to
  \\k\\ derivatives. The method is never routed to quadrature at any
  order.

## Defaults of each function

When `method` is `NULL`, the function returns the closed form of the
family where there is one. The route used otherwise is:

|  |  |
|----|----|
| [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md) | [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md) |
| [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md) and [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md) with `expected = TRUE` | [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md) (`identity = "hessian"`) |
| [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md) and [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md) | [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md) (order 2) |
| [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md), [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md), [`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md) and [`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md) | none: a family without a closed form signals an error |

These defaults reproduce the values that the functions returned before
`method` was introduced.

## Derivatives of the expected information

The functions that return the derivatives of the expected information
with respect to the parameters have the closed form of the family as
their only default. With a method, the \\m\\th derivative of the matrix
\\E\_{ab}(\theta)\\ that the method defines is returned.

For
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
the matrix is the expectation \\E\_{ab}(\theta) =
\mathbb{E}\_\theta\[f\_{ab}\]\\ of the integrand \\f\_{ab}\\ of the
identity (\\\ell\_{ab}\\ for the Hessian identity and \\-\ell_a\ell_b\\
for the score identity), and the derivative moves the density as well as
the integrand,
\$\$\frac{\partial}{\partial\theta_c}\mathbb{E}\_\theta\[f\] =
\mathbb{E}\_\theta\left\[\frac{\partial f}{\partial\theta_c} +
f\\\ell_c\right\],\$\$ repeated for higher orders. The expectation of
the right-hand side is taken by the method.

For
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
the matrix is the sum over the observations of \\-\ell_a\ell_b\\ at the
data, which is the matrix whose inverse the scoring step uses, and the
density does not move: the exact derivative of that matrix with respect
to \\\theta_c\\ is, per observation, \$\$-\left(\ell\_{ac}\\\ell_b +
\ell_a\\\ell\_{bc}\right),\$\$ and the second derivative is
\$\$-\left(\ell\_{acd}\\\ell_b + \ell\_{ac}\\\ell\_{bd} +
\ell\_{ad}\\\ell\_{bc} + \ell_a\\\ell\_{bcd}\right).\$\$ These are not
estimates of \\\partial \mathbb{E}\[\ell\_{ab}\]/\partial \theta_c\\,
which contains the moved density; they are the derivatives of the matrix
that
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
returns, and a criterion built on that matrix needs these. The
derivatives of order 3 and 4 follow the same rule by Leibniz's formula.

A derivative of order \\m\\ of the matrix needs the observed derivatives
of the log-density up to order \\m+1\\ for the score identity, which is
always the case for \\m \le 4\\, and up to order \\m+2\\ for the Hessian
identity, which is available for \\m \le 3\\; for \\m = 4\\ the Hessian
identity signals an error that names `identity = "score"`.

## A kink in the log-density

Where the log-density has a kink in a parameter, the observed second
derivative has a point mass at the kink, and an average of its pointwise
values does not contain that mass. For the Laplace location,
\\\ell\_{\mu\mu}\\ is zero almost everywhere, so its pointwise average
is zero, while the information is \\1/\sigma^2\\. The information of
such a family is the variance of the score, which the score identity
evaluates. For a family that declares a non-smooth parameter (see
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)),
the Hessian identity at order 2, whether by
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
or by
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md),
signals an error that names `identity = "score"`.
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
which uses the score identity, and `identity = "score"` remain valid.
When the family declares where its kink is, through
[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md),
the quadrature of
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
places a knot there.

## See also

[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
for the constructors;
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
for the classes;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
the generics that take `method`;
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
for the derivatives of the expected information;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
where it is set for a fit;
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
and
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md),
the objects that select the information of a fit, and
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
for their comparison;
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md),
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
and
[`expected_by_mc()`](https://statmodels7.github.io/distributions7/reference/expected_by_mc.md)
for the implementations.

## Examples

``` r
# A family defined by its density alone has no closed-form information.
Dens <- S7::new_class("Dens", parent = continuous_distrib)
S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(y, theta[[1]], theta[[2]], log = log)
}
d <- Dens(distrib_name = "dens", dimension = "univariate",
          bounds = c(-Inf, Inf), params = c("mu", "sigma"),
          params_interpretation = c(mu = "location", sigma = "scale"),
          n_params = 2L,
          params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
          link_params = list(mu = linkfunctions7::identity_link(),
                             sigma = linkfunctions7::log_link()),
          params_smooth = c(mu = TRUE, sigma = TRUE))
th <- list(mu = 0, sigma = 2)
y <- c(-1, 0.5, 2)

# The sigma-sigma entry of the expected Hessian is -2 / sigma^2 = -0.5.
# The two quadratures agree; the empirical outer product of the scores is
# the value at one observation, and Monte Carlo carries sampling error.
methods <- list(
  opg = expect_opg(),
  integrate_hessian = expect_integrate(identity = "hessian"),
  integrate_score = expect_integrate(identity = "score"),
  mc_hessian = expect_mc(nsim = 2000, identity = "hessian"),
  mc_score = expect_mc(nsim = 2000, identity = "score")
)
set.seed(3)
vapply(methods, function(m) {
  distrib_expected_hessian(d, y, th, method = m)$sigma_sigma[1]
}, numeric(1))
#>               opg integrate_hessian   integrate_score        mc_hessian 
#>        -0.1406250        -0.5000000        -0.5000000        -0.4575573 
#>          mc_score 
#>        -0.4840274 

# An object is used also by a family that has a closed form.
g <- gaussian1_distrib()
distrib_expected_hessian(g, 0, list(mu = 0, sigma = 2))$sigma_sigma
#> [1] -0.5
distrib_expected_hessian(g, 0, list(mu = 0, sigma = 2),
                         method = expect_integrate())$sigma_sigma
#> [1] -0.5
```
