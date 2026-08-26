# Skew t Score

Computes the four first derivatives of the log-density. Three are closed
form: with \\D = A + QB\\ in the notation of
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md),
\$\$\dfrac{\partial \ell}{\partial \mu} = -\dfrac{D}{\sigma}, \qquad
\dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1 + zD}{\sigma}, \qquad
\dfrac{\partial \ell}{\partial \alpha} = Q z c.\$\$

The fourth is not. \\\partial\log T\_{\nu+1}(w)/\partial\nu\\ is a
derivative of a Student \\t\\ distribution function with respect to its
degrees of freedom, which has no elementary expression, the same
obstruction the gamma and beta distribution functions meet in their
shape. That one component is a single central difference of the
**log-density**, taken with
[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
at the step of
[`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md).

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
  numeric vector of length 1 or of the length of `y`. `sigma` and `nu`
  must be strictly positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of four numeric vectors, `mu`, `sigma`, `alpha` and `nu`,
each of the length of the recycled inputs.

## Accuracy

Measured at \\\mu = 0\\, \\\sigma = 1\\, \\\alpha = 3\\, \\\nu = 6\\ on
four observations, the summed score agrees with
[`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html) on the
log-likelihood to \\6\times10^{-12}\\ in \\\mu\\, \\2\times10^{-11}\\ in
\\\sigma\\, \\3\times10^{-12}\\ in \\\alpha\\ and \\5\times10^{-11}\\ in
\\\nu\\. The \\\nu\\ component is therefore the loosest of the four, and
it is what sets the accuracy a fit of this family can stop at; see
[`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md).

## Notation

\\z = (y-\mu)/\sigma\\, \\c = \sqrt{(\nu+1)/(\nu+z^2)}\\, \\w = \alpha z
c\\, and \\A\\, \\B\\, \\Q\\ are as
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
defines them.

## See also

[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
for the scalar functions,
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
for the second derivatives,
[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
for the stencil, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
g <- distrib_gradient(d, y, th)

# The location component written out in skewt_pieces()' terms.
p <- distributions7:::skewt_pieces(y, 0, 1, 3, 6)
all.equal(g$mu, -(p$a + p$q * p$b))
#> [1] TRUE

# All four against numerical differentiation of the log-likelihood.
f <- function(v) sum(distrib_pdf(d, y, list(mu = v[1], sigma = v[2],
                                            alpha = v[3], nu = v[4]),
                                 log = TRUE))
rbind(analytic = vapply(g, sum, 0),
      numeric = numDeriv::grad(f, c(0, 1, 3, 6)))
#>                 mu    sigma     alpha         nu
#> analytic -6.892263 5.904851 -2.056567 -0.3632304
#> numeric  -6.892263 5.904851 -2.056567 -0.3632304

# At shape zero the score in alpha is not zero: the tilting factor is at
# its inflection, so alpha is still identified.
distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 0, nu = 6))$alpha
#> [1] -1.0638843 -0.2476525  0.3283218  1.3259405
```
