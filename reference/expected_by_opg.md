# The Information as the Outer Product of the Observed Scores

Returns \\\ell^{(ij)} \approx -\ell^{(i)}\ell^{(j)}\\ evaluated at each
observation, the per-observation outer product of the score. It is the
BHHH estimator of the information, and unlike every other route here it
takes no expectation: it integrates nothing and sums nothing over the
support, so it costs one call to
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
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md);
each has the length of `y`.

## Details

This route differs from the identity from which it comes. The second
Bartlett identity states that \\\mathbb{E}\[\ell^{(ij)}\] =
-\mathbb{E}\[\ell^{(i)}\ell^{(j)}\]\\, and
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
evaluates the right-hand side as written, which for a discrete family is
a sum over the whole support and for a continuous one a quadrature. This
function drops the expectation and reads the integrand at the
observation in hand. The two agree in expectation and differ in cost:
for a discrete family the identity evaluates the derivative kernel over
the whole support for every observation, whereas this route evaluates it
once per observation.

The estimate has more variance than the exact expectation and remains
consistent. \\-\ell^{(i)}\ell^{(j)}\\ is an unbiased estimate of the
corresponding component of \\\mathbb{E}\[\ell^{(ij)}\]\\ at that
observation, so the sum over observations, which a scoring step
aggregates into \\X^\top W X\\, estimates the information consistently.
Component by component it is a rank-one matrix and correlates poorly
with the exact expectation; summed over the observations it approximates
the exact information, and a fit driven by it reaches the same maximum.

It is the right default for a scoring step because it is a sum of outer
products and is therefore positive semidefinite by construction, whereas
the observed Hessian need not be positive semidefinite away from the
optimum, and because the fixed point is unchanged, since the score is
exact and any positive definite matrix takes a scoring iteration to the
same stationary point. The route is not suited to a reported standard
error, which is read off the observed information instead; see
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
and its [`vcov()`](https://rdrr.io/r/stats/vcov.html) method.

The function applies at order 2 only. The outer product of scores is the
second-order identity and has no counterpart above it, so
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
routes a higher order to
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md).

## See also

[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
for the identity that this function estimates,
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the catalog of routes, and
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the predicate that tests whether a family needs any of them.
