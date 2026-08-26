# Strategies for Expected Derivatives

When a distribution does not supply a closed-form expected derivative,
the expectation has to be approximated. The generics
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
therefore accept an `approx` argument selecting *how* the expectation is
taken. The argument is **ignored** when the distribution provides an
analytical method (which is always preferred).

## Value

Nothing. This page documents the `approx` argument shared by the three
generics named above; the value returned is theirs.

## Details

Let \\\ell\\ be the log-density and \\\ell\_{i}, \ell\_{ij}, \dots\\ its
derivatives with respect to the parameters. The three strategies are:

- `"bartlett"` (also spelled `"opg"`):

  Uses the Bartlett identity of the corresponding order, which expresses
  the expected derivative through expectations of *products of
  lower-order* derivatives. Writing the identity as a sum over set
  partitions of the index set, \$\$\sum\_{\pi}
  \mathbb{E}\left\[\prod\_{B \in \pi} \ell_B\right\] = 0,\$\$ the target
  term (the single-block partition) is obtained from all the others. At
  order 2 this is exactly the familiar **outer product of gradients**,
  \\\mathbb{E}\[\ell\_{ij}\] = -\mathbb{E}\[\ell_i \ell_j\]\\; at order
  3, \\\mathbb{E}\[\ell\_{ijk}\] =
  -\left(\mathbb{E}\[\ell\_{ij}\ell_k\] +
  \mathbb{E}\[\ell\_{ik}\ell_j\] + \mathbb{E}\[\ell\_{jk}\ell_i\] +
  \mathbb{E}\[\ell_i\ell_j\ell_k\]\right)\\, and so on. *Fastest at
  order 2* (only the score is needed) and the only variant that stays
  valid when the log-likelihood is not differentiable in a parameter,
  where \\\mathbb{E}\[H\]\\ degenerates but the score variance still
  gives the information (see
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)).
  Deterministic. At higher orders it needs several integrals, so it is
  usually slower than `"integrate"`. It never requires the top-order
  derivative itself, which is useful where only the lower orders are
  available in closed form.

- `"integrate"`:

  Integrates the observed derivative of that order directly against the
  density (numerical quadrature for continuous distributions, series
  summation for discrete ones). Deterministic and normally the most
  accurate when the observed derivative is available in closed form.
  Estimates \\\mathbb{E}\[\partial^k \ell\]\\ literally, which for a
  non-regular model is *not* the information.

- `"mc"`:

  Simulates `nsim` observations from the distribution and averages the
  observed derivative over them. The simplest and most robust option
  when quadrature struggles (heavy tails, awkward supports), but
  stochastic: the error decreases only as \\1/\sqrt{n\_{sim}}\\, so
  results are not exactly reproducible unless the seed is fixed.
  Estimates the same quantity as `"integrate"`. **Note:** the cost of
  this option is the cost of simulating from the distribution. A native
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
  is best, and failing that the default RNG falls back to
  [`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md),
  which only needs the density. It becomes expensive only for a
  distribution that supplies a quantile function slow enough to make
  inverse transform sampling the bottleneck; prefer `"bartlett"` or
  `"integrate"` in that case.

**Defaults.** `distrib_expected_hessian` defaults to `"bartlett"`,
because at order 2 it is both the cheapest (only first derivatives) and
the most broadly valid. `distrib_deriv3` and `distrib_deriv4` default to
`"integrate"`, since at those orders direct integration of the available
derivative is usually cheaper and more accurate.

**What the kink costs, measured.** On a Laplace carrying a density, a
score and a Hessian but no expected method, at \\\sigma = 1\\ over 200
observations: `"bartlett"` returns \\-200\\, which is \\-n/\sigma^2\\
and agrees with the shipped family's closed form to the digit, while
`"integrate"` and `"mc"` both return **exactly 0**. Neither is wrong
about what it computes. The observed \\\ell\_{\mu\mu}\\ really is zero
almost everywhere, so its expectation is zero; what fails is the
identification of that expectation with \\-\mathcal{I}(\theta)\\, which
is the second Bartlett identity. Only the score-based route survives,
and the information of a non-regular family is *defined* as the variance
of the score.

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
# A family with no closed-form expected information reads 'approx'. The two
# deterministic strategies agree; Monte Carlo is the same quantity with
# sampling error on it.
sn <- skewnormal1_distrib()
th <- list(mu = 0, sigma = 1, alpha = 3)
set.seed(2)
y <- distrib_rng(sn, 40, th)

vapply(c("bartlett", "integrate"), function(a) {
  sum(distrib_expected_hessian(sn, y, th, approx = a)$alpha_alpha)
}, numeric(1))
#>  bartlett integrate 
#> -1.085353 -1.085353 

set.seed(3)
sum(distrib_expected_hessian(sn, y, th, approx = "mc", nsim = 500)$alpha_alpha)
#> [1] -1.064756

# A family that writes its expected information out ignores the argument,
# and fit_distrib() rejects one given there rather than dropping it.
g <- gaussian1_distrib()
identical(distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1)),
          distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1),
                                   approx = "mc"))
#> [1] TRUE
```
