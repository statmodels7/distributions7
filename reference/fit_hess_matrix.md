# Summed Hessian on the Link Scale, as a Matrix

Assembles the package's named list of second-derivative components into
the symmetric \\p \times p\\ matrix that an optimizer requires, summed
over observations and on the link scale. Components are stored one per
unordered index pair, so this fills both triangles from the one value.

With `expected = TRUE` it assembles the expected Hessian instead, which
turns Newton's method into Fisher scoring and makes a fit possible on a
family whose observed curvature is unusable. On a Laplace at \\\sigma =
1\\ with 400 draws the observed matrix is \\\bigl(\begin{smallmatrix} 0
& 2\\ 2 & -392\end{smallmatrix}\bigr)\\:
\\\partial^2\ell/\partial\mu^2\\ is zero almost everywhere, so the
determinant is \\-(\sum_i \mathrm{sign}(y_i - \mu))^2/\sigma^2\\, which
is negative unless the signs balance exactly, and the matrix is
indefinite. The expected information at the same point is \\-400 I\\.

## Usage

``` r
fit_hess_matrix(distrib, y, theta, expected, method = NULL, threads = 1L)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations, or the response matrix of a
  multivariate family.

- theta:

  A named list of parameters on the parameter scale, aligned to
  `distrib@params`.

- expected:

  Logical of length 1, with no default. `TRUE` assembles the expected
  Hessian and `FALSE` the observed one.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
  that selects how the expected information is obtained: `NULL` takes
  the closed form of the family where it has one and
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  otherwise, and an object is used for every family; see
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).
  It is read only when `expected` is `TRUE`.

- threads:

  How many threads the family's compiled kernels may use, as an integer
  count. Defaults to `1L`. The matrix does not depend on the count.

## Value

A symmetric numeric matrix of dimension `length(distrib@params)`, whose
dimnames are both `distrib@params`.

## See also

[`fit_score()`](https://statmodels7.github.io/distributions7/reference/fit_score.md)
for the first-order counterpart;
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
which supply the components;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
where `method` is set.
