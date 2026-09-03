# The Information as the Outer Product of the Observed Scores

Returns \\\ell^{(ij)} \approx -\ell^{(i)}\ell^{(j)}\\ evaluated at each
observation, the per-observation outer product of the score. It is the
BHHH estimator of the information, and unlike every other route here it
takes no expectation: nothing is integrated and nothing is summed over
the support, so it costs one call to
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md).

## Usage

``` r
expected_by_opg(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list of expected Hessian component vectors, named by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each of the length of `y`.

## Details

The name is worth keeping apart from the identity it comes from. The
second Bartlett identity states that \\\mathbb{E}\[\ell^{(ij)}\] =
-\mathbb{E}\[\ell^{(i)}\ell^{(j)}\]\\, and
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
evaluates the right-hand side as written, which for a discrete family is
a sum over the whole support and for a continuous one a quadrature. This
function drops the expectation and reads the integrand at the
observation in hand. The two agree in expectation and not in cost:
measured on
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
with \\\mu\\ varying by observation, one call to the identity at \\n =
1000\\ evaluated 1.67 million rows of the family's derivative kernel and
took 18.4 seconds, where this route evaluates \\n\\ rows.

What is given up is variance, not consistency. \\-\ell^{(i)}\ell^{(j)}\\
is an unbiased estimate of the corresponding component of
\\\mathbb{E}\[\ell^{(ij)}\]\\ at that observation, so the sum over
observations, which is what a scoring step aggregates into \\X^\top W
X\\, estimates the information consistently. Component by component it
is a rank-one matrix and correlates poorly with the exact expectation;
summed it agreed with the exact route to within 10 per cent on the same
measurement, and a fit driven by it reached the same maximum.

Two properties make it the right default for a scoring step. It is
positive semidefinite by construction, being a sum of outer products,
where the observed Hessian need not be away from the optimum; and the
fixed point is unchanged, because the score is exact and any positive
definite matrix takes a scoring iteration to the same stationary point.
What it is not suited to is a reported standard error, which is read off
the observed information instead – see
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
and its [`vcov()`](https://rdrr.io/r/stats/vcov.html) method.

Order 2 only. The outer product of scores is the second-order identity
and has no counterpart above it, so
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
routes a higher order to
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md).

## See also

[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
for the identity this estimates,
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the catalog of routes, and
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the predicate that says whether a family needs any of them.
