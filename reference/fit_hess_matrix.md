# Summed Hessian on the Link Scale, as a Matrix

Assembles the package's named list of second-derivative components into
the symmetric \\p \times p\\ matrix an optimizer wants, summed over
observations and on the link scale. Components are stored one per
unordered index pair, so this fills both triangles from the one value.

With `expected = TRUE` it assembles the expected Hessian instead, which
is what turns Newton's method into Fisher scoring and what makes a fit
possible on a family whose observed curvature is unusable. On a Laplace
at \\\sigma = 1\\ with 400 draws the observed matrix is
\\\bigl(\begin{smallmatrix} 0 & 2\\ 2 & -392\end{smallmatrix}\bigr)\\:
\\\partial^2\ell/\partial\mu^2\\ is zero almost everywhere, so the
determinant is \\-(\sum_i \mathrm{sign}(y_i - \mu))^2/\sigma^2\\, which
is negative unless the signs balance exactly, and the matrix is
**indefinite**. The expected information at the same point is \\-400
I\\.

## Usage

``` r
fit_hess_matrix(
  distrib,
  y,
  theta,
  expected,
  approx = "bartlett",
  nsim = 10000,
  threads = 1L
)
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

- approx:

  How the expectation is approximated for a family with no closed form
  for it: `"bartlett"` (the default), `"integrate"` or `"mc"`. Read only
  when `expected` is `TRUE` and the family has no exact expression;
  ignored otherwise.

- nsim:

  Monte Carlo sample size, a single positive integer, read only when
  `approx = "mc"` is in force. Defaults to 10000.

- threads:

  How many threads the family's compiled kernels may use, as an integer
  count. Defaults to `1L`. The matrix does not depend on the count.

## Value

A symmetric numeric matrix of dimension `length(distrib@params)`, with
both dimnames set to `distrib@params`.

## See also

[`fit_score()`](https://statmodels7.github.io/distributions7/reference/fit_score.md)
for the first-order counterpart;
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
which supply the components;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
where `approx` and `nsim` are set.
