# The Empirical Outer Product of the Scores

Returns the score identity of
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
evaluated at each observation, with no expectation. At order 2 this is
\\-\ell^{(i)}\ell^{(j)}\\, the per-observation outer product of the
score (the BHHH estimator of the information); at orders 3 and 4 it is
the sum of products \\-\sum\_{\pi,\\\|\pi\| \ge 2} \prod\_{B \in \pi}
\ell_B\\. Unlike every other route it takes no expectation: it
integrates nothing and sums nothing over the support, so it costs the
derivatives of orders 1 to `order - 1` at the observations.

## Usage

``` r
expected_by_opg(distrib, y, theta, order = 2L)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4. The default is 2.

## Value

A named list of expected derivative component vectors, named as
[`observed_deriv()`](https://statmodels7.github.io/distributions7/reference/observed_deriv.md)
names that order; each has the length of `y`.

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
consistent: each term is an unbiased estimate of the corresponding
component of the expected derivative at that observation, so the sum
over observations estimates it consistently. Component by component at
order 2 it is a rank-one matrix and correlates poorly with the exact
expectation; summed over the observations it approximates the exact
information, and a fit driven by it reaches the same maximum. At orders
3 and 4 the terms are products of up to four derivatives, and the
variance is larger still.

It is the right default for a scoring step because it is a sum of outer
products and is therefore positive semidefinite by construction, whereas
the observed Hessian need not be positive semidefinite away from the
optimum, and because the fixed point is unchanged, since the score is
exact and any positive definite matrix takes a scoring iteration to the
same stationary point. The route is not suited to a reported standard
error, which is read off the observed information instead; see
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
and its [`vcov()`](https://rdrr.io/r/stats/vcov.html) method.

## See also

[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
for the identity that this function estimates,
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the catalog of routes, and
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the predicate that tests whether a family needs any of them.
