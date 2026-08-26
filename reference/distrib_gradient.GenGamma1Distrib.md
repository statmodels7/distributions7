# Generalized Gamma Score

Computes the three first derivatives of the log-density in closed form,
in a compiled kernel. With \\w = (y/a)^{p}\\, \\L = \log(y/a)\\ and \\k
= d/p\\, \$\$\dfrac{\partial\ell}{\partial a} = \dfrac{pw - d}{a},
\qquad \dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
\dfrac{\partial\ell}{\partial p} = \dfrac{1}{p} +
\dfrac{d\\\psi(k)}{p^{2}} - wL.\$\$ Only the digamma function appears
beyond elementary operations, and it appears at one argument, \\k =
d/p\\.

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

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of three numeric vectors, `a`, `d` and `p`, each of the
length of the recycled inputs.

## Notation

\\\ell\\ is the log-density of one observation, \\\psi\\ the digamma
function [`digamma()`](https://rdrr.io/r/base/Special.html), and \\a\\,
\\d\\, \\p\\ the scale and the two shapes.

## See also

[`distrib_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GenGamma1Distrib.md)
for the closed-form information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
y <- c(0.5, 1.5, 4)
th <- list(a = 2, d = 3, p = 1.5)
g <- distrib_gradient(d, y, th)

# The scale component written out.
all.equal(g$a, (1.5 * (y / 2)^1.5 - 3) / 2)
#> [1] TRUE

# All three against numerical differentiation of the log-likelihood.
f <- function(v) sum(distrib_pdf(d, y, list(a = v[1], d = v[2], p = v[3]),
                                 log = TRUE))
rbind(analytic = vapply(g, sum, 0),
      numeric = numDeriv::grad(f, c(2, 3, 1.5)))
#>                 a         d        p
#> analytic -1.79779 -1.826398 2.090763
#> numeric  -1.79779 -1.826398 2.090763

# The score sums to nearly zero at the maximum likelihood estimate.
set.seed(52)
x <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, x)
vapply(distrib_gradient(d, x, as.list(coef(fit))), sum, 0) / 2000
#>             a             d             p 
#>  2.063203e-09 -2.914640e-09 -1.328359e-08 
```
