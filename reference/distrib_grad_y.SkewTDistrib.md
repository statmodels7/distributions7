# Skew t Response Derivative

Computes \\\partial\ell/\partial y = D/\sigma\\, with \\D = A + QB\\ in
the notation of
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md).
It is minus the derivative in \\\mu\\, exactly, because the response and
the location enter the density only through their difference; the
identity holds for every location family.

It is closed form at every parameter value, \\\nu\\ included: nothing
here differentiates the degrees of freedom.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
  numeric vector of length 1 or of the length of `y`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the length of the recycled inputs.

## Notation

\\\ell\\ is the log-density of one observation and \\A\\, \\B\\, \\Q\\
are as
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
defines them.

## See also

[`distrib_hess_y.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewTDistrib.md)
for the second derivative,
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md)
for the parameter derivatives, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)

# It is minus the score in the location, exactly.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Against a central difference of the log-density.
eps <- 1e-6
rbind(analytic = distrib_grad_y(d, y, th),
      numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
                 distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#>              [,1]    [,2]      [,3]      [,4]
#> analytic 3.787052 4.38216 0.1341366 -1.411086
#> numeric  3.787052 4.38216 0.1341366 -1.411086

# The score redescends in the heavy tail, as a Student t's does.
distrib_grad_y(d, c(2, 8, 32, 128), th)
#> [1] -1.39870409 -0.79999568 -0.21747567 -0.05466748
```
