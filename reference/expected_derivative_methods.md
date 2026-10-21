# Strategies for Expected Derivatives

When a distribution does not supply a closed-form expected derivative,
the expectation has to be approximated. The generics
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
therefore accept an `approx` argument selecting *how* the expectation is
taken. The argument is ignored when the distribution provides an
analytical method (which is always preferred).

## Value

Nothing. This page documents the `approx` argument shared by the three
generics named above; the value returned is theirs.

## Details

Let \\\ell\\ be the log-density and \\\ell\_{i}, \ell\_{ij}, \dots\\ its
derivatives with respect to the parameters. The strategies are:

- `"opg"`:

  At order 2, reads \\-\ell_i \ell_j\\ at each observation and takes no
  expectation, at the cost of one score evaluation; see
  [`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md).
  Above order 2 it is routed to `"bartlett"`.

- `"bartlett"`:

  Uses the Bartlett identity of the corresponding order, which expresses
  the expected derivative through expectations of *products of
  lower-order* derivatives. Writing the identity as a sum over set
  partitions of the index set, \$\$\sum\_{\pi}
  \mathbb{E}\left\[\prod\_{B \in \pi} \ell_B\right\] = 0,\$\$ the target
  term (the single-block partition) is obtained from all the others. At
  order 2 this is exactly the familiar outer product of gradients,
  \\\mathbb{E}\[\ell\_{ij}\] = -\mathbb{E}\[\ell_i \ell_j\]\\; at order
  3, \\\mathbb{E}\[\ell\_{ijk}\] =
  -\left(\mathbb{E}\[\ell\_{ij}\ell_k\] +
  \mathbb{E}\[\ell\_{ik}\ell_j\] + \mathbb{E}\[\ell\_{jk}\ell_i\] +
  \mathbb{E}\[\ell_i\ell_j\ell_k\]\right)\\, and so on. At order 2 it
  needs only the score, whose product it integrates or sums, so it is
  slower than `"opg"`. It stays valid when the log-likelihood is not
  differentiable in a parameter, where \\\mathbb{E}\[H\]\\ degenerates
  but the score variance still gives the information (see
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)).
  It is deterministic. At higher orders it needs several integrals, so
  it is usually slower than `"integrate"`. It never requires the
  top-order derivative itself, which is useful where only the lower
  orders are available in closed form.

- `"integrate"`:

  Integrates the observed derivative of that order directly against the
  density (numerical quadrature for continuous distributions, series
  summation for discrete ones). It is deterministic and normally the
  most accurate when the observed derivative is available in closed
  form. It estimates \\\mathbb{E}\[\partial^k \ell\]\\ literally, with
  one exception: at order 2, for a family with a non-smooth parameter
  (see
  [`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)),
  it evaluates \\-\mathbb{E}\[\ell_i \ell_j\]\\ as `"bartlett"` does,
  because the observed \\\ell\_{ij}\\ then lacks the point mass at the
  kink.

- `"mc"`:

  Simulates `nsim` observations from the distribution and averages the
  observed derivative over them. The simplest and most robust option
  when quadrature performs poorly (heavy tails, awkward supports), but
  it is stochastic: the error decreases only as \\1/\sqrt{n\_{sim}}\\,
  so results are not exactly reproducible unless the seed is fixed.
  Estimates the same quantity as `"integrate"`. The cost of this option
  is the cost of simulating from the distribution. A native
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
  is best; failing that, the default generator uses inverse transform
  when the family has a quantile method of its own and
  [`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md),
  which only needs the density, otherwise. It becomes expensive only for
  a distribution that supplies a quantile function slow enough to make
  inverse transform sampling the bottleneck; prefer `"bartlett"` or
  `"integrate"` in that case.

The generic `distrib_expected_hessian` defaults to `"opg"`, which needs
one score evaluation per observation and no expectation.
`distrib_deriv3` and `distrib_deriv4` default to `"integrate"`, since at
those orders direct integration of the available derivative is usually
cheaper and more accurate.

Where the log-density has a kink in a parameter, the observed second
derivative has a point mass at the kink, and an average of its pointwise
values does not contain that mass. For the Laplace location,
\\\ell\_{\mu\mu}\\ is zero almost everywhere, so its pointwise average
is zero, while the information is \\1/\sigma^2\\. The information of
such a family is the variance of the score, and the second Bartlett
identity evaluates that variance. `"bartlett"`, `"opg"` and, at order 2,
`"integrate"` therefore return \\-1/\sigma^2\\ per observation; `"mc"`
averages the observed \\\ell\_{\mu\mu}\\ and returns 0. When the family
declares where its kink is, through
[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md),
the quadrature of
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
places a knot there.

## See also

[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
the three generics that take `approx`;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
where it is set for a fit;
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md),
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
and
[`expected_by_mc()`](https://statmodels7.github.io/distributions7/reference/expected_by_mc.md)
for the three implementations.

## Examples

``` r
# A family defined by its density alone reads 'approx'. The two
# deterministic strategies agree on the exact -2 / sigma^2; "opg" reads one
# observation's score product, and Monte Carlo carries sampling error.
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

vapply(c("opg", "bartlett", "integrate"), function(a) {
  distrib_expected_hessian(d, y, th, approx = a)$sigma_sigma[1]
}, numeric(1))
#>       opg  bartlett integrate 
#> -0.140625 -0.500000 -0.500000 

set.seed(3)
distrib_expected_hessian(d, y, th, approx = "mc", nsim = 2000)$sigma_sigma[1]
#> [1] -0.4575573

# A family that writes its expected information out ignores the argument,
# and fit_distrib() rejects one given there rather than dropping it.
g <- gaussian1_distrib()
identical(distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1)),
          distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1),
                                   approx = "mc"))
#> [1] TRUE
```
