# Generalized Pareto Observed Hessian

Computes the three second derivatives of the log-density in closed form,
in a compiled kernel. In the notation of
[`distrib_gradient.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GPDDistrib.md)
the expressions stay short because \\t - \xi z = 1\\, which makes
\\\partial u/\partial\sigma\\ equal to \\-z/(\sigma t^2)\\.

The pure-\\\xi\\ component carries the same removable singularity the
score does, one order worse, and goes through the same analytic function
\\\Lambda(u) = \log(1+u)/u\\. The two terms that individually diverge as
\\\xi \to 0\\ cancel, and the limit is an ordinary point of the formula.

This is the **observed** curvature at the data. The expected one is
closed form as well, but only above \\\xi = -1/2\\; see
[`distrib_expected_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GPDDistrib.md).

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of three numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `sigma_sigma`, `xi_xi`, `sigma_xi`.

## Notation

\\z = y/\sigma\\, \\t = 1 + \xi z\\, \\u = z/t\\, and \\\ell\\ is the
log-density of one observation.

## See also

[`distrib_gradient.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GPDDistrib.md)
for the order below,
[`distrib_deriv3.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
for the order above,
[`distrib_expected_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GPDDistrib.md)
for the expectation, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
y <- c(0.2, 1, 4)
th <- list(sigma = 1.5, xi = 0.3)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "sigma_sigma" "sigma_xi"    "xi_xi"      

# Against a central difference of the score.
eps <- 1e-5
rbind(analytic = h$sigma_xi,
      numeric = (distrib_gradient(d, y, list(sigma = 1.5,
                                             xi = 0.3 + eps))$sigma -
                 distrib_gradient(d, y, list(sigma = 1.5,
                                             xi = 0.3 - eps))$sigma) /
                (2 * eps))
#>                [,1]      [,2]       [,3]
#> analytic 0.07122507 0.1028807 -0.9144947
#> numeric  0.07122507 0.1028807 -0.9144947

# The pure-shape component passes through zero without a branch.
vapply(c(1e-2, 1e-6, 0, -1e-6),
       function(x) distrib_hessian(d, 1, list(sigma = 1, xi = x))$xi_xi, 0)
#> [1] 0.3283927 0.3333328 0.3333333 0.3333338
```
