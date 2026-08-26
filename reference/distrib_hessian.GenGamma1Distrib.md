# Generalized Gamma Observed Hessian

Computes the six second derivatives of the log-density in closed form,
in a compiled kernel, by differentiating the expressions of
[`distrib_gradient.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md)
again.

One component is worth naming: the mixed \\a\\-\\d\\ entry is \\-1/a\\,
free of the data. The scale and the first shape meet in the log-density
only through the term \\-d\log a\\, which is bilinear, so its second
derivative carries no observation.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- y:

  A numeric vector of observations, strictly positive.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`. All three must be strictly
  positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of six numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `a_a`, `d_d`, `p_p`, `a_d`, `a_p`, `d_p`.

## Notation

\\w = (y/a)^p\\, \\L = \log(y/a)\\, \\k = d/p\\, and \\\psi\\, \\\psi'\\
are the digamma and trigamma functions.

## See also

[`distrib_gradient.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md)
for the order below,
[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
for the order above,
[`distrib_expected_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GenGamma1Distrib.md)
for its expectation, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
y <- c(0.5, 1.5, 4)
th <- list(a = 2, d = 3, p = 1.5)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "a_a" "a_d" "a_p" "d_d" "d_p" "p_p"

# The scale-shape entry is -1/a at every observation.
h$a_d
#> [1] -0.5 -0.5 -0.5

# Against a central difference of the score.
eps <- 1e-5
rbind(analytic = h$d_p,
      numeric = (distrib_gradient(d, y, list(a = 2, d = 3,
                                             p = 1.5 + eps))$d -
                 distrib_gradient(d, y, list(a = 2, d = 3,
                                             p = 1.5 - eps))$d) / (2 * eps))
#>               [,1]      [,2]      [,3]
#> analytic 0.7611789 0.7611789 0.7611789
#> numeric  0.7611789 0.7611789 0.7611789
```
