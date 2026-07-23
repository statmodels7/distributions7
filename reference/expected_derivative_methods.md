# Strategies for Expected Derivatives

When a distribution does not supply a closed-form expected derivative,
the expectation has to be approximated. The generics
[`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
therefore accept an `approx` argument selecting *how* the expectation is
taken. The argument is **ignored** when the distribution provides an
analytical method (which is always preferred).

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
  [`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)).
  Deterministic. At higher orders it needs several integrals, so it is
  usually slower than `"integrate"`, but it never requires the top-order
  derivative itself – useful when only the lower orders are available in
  closed form.

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
  [`distrib_rng`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
  is best, and failing that the default RNG falls back to
  [`rng_grou`](https://statmodels7.github.io/distributions7/reference/rng_grou.md),
  which only needs the density. It becomes expensive only for a
  distribution that supplies a quantile function slow enough to make
  inverse transform sampling the bottleneck; prefer `"bartlett"` or
  `"integrate"` in that case.

**Defaults.** `distrib_expected_hessian` defaults to `"bartlett"`,
because at order 2 it is both the cheapest (only first derivatives) and
the most broadly valid. `distrib_deriv3` and `distrib_deriv4` default to
`"integrate"`, since at those orders direct integration of the available
derivative is usually cheaper and more accurate.

## See also

[`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
