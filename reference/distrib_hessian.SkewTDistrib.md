# Skew t Observed Hessian

Computes the ten second derivatives of the log-density. The block in
\\(\mu, \sigma, \alpha)\\ is closed form; every component involving
\\\nu\\ comes from one stencil applied to an analytic quantity, for the
reason
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md)
gives.

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

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body, so this method always returns the
  parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of ten numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `mu_mu`, `sigma_sigma`, `alpha_alpha`, `nu_nu`, then `mu_sigma`,
`mu_alpha`, `mu_nu`, `sigma_alpha`, `sigma_nu`, `alpha_nu`.

## The closed-form block

With \\D = A + QB\\ and \\D' = A' + Q'B^2 + QB'\\ in the notation of
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md),
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{D'}{\sigma^2},
\qquad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial \sigma} =
\dfrac{D + zD'}{\sigma^2}, \qquad \dfrac{\partial^2 \ell}{\partial
\sigma^2} = \dfrac{1 + 2zD + z^2 D'}{\sigma^2},\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \alpha^2} = Q' z^2 c^2, \qquad
\dfrac{\partial^2 \ell}{\partial \mu \\ \partial \alpha} = -\dfrac{Q' B
z c + Q E}{\sigma}, \qquad \dfrac{\partial^2 \ell}{\partial \sigma \\
\partial \alpha} = -\dfrac{z(Q' B z c + Q E)}{\sigma}.\$\$

## The four components in the degrees of freedom

\\\partial^2\ell/\partial\nu^2\\ is one five-point second difference of
the log-density,
[`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md).
The three mixed ones step the **closed-form score** in \\\nu\\ with
[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md),
so only one difference is taken and it is taken of an analytic quantity.
Stepping the log-density in both variables instead would be a difference
of a difference in \\\nu\\, which this family never does.

Measured at \\\mu = 0\\, \\\sigma = 1\\, \\\alpha = 3\\, \\\nu = 6\\,
the summed `nu_nu` agrees with
[`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html) on
the log-likelihood to \\2\times10^{-9}\\ relative and `mu_nu` to the
printed digit.

## Notation

\\z = (y-\mu)/\sigma\\, \\c = \sqrt{(\nu+1)/(\nu+z^2)}\\, and \\A\\,
\\B\\, \\E\\, \\Q\\ are as
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
defines them.

## See also

[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md)
for the order below,
[`distrib_deriv3.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md)
for the order above,
[`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md)
for the stencil, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
h <- distrib_hessian(d, y, th)
names(h)
#>  [1] "mu_mu"       "sigma_sigma" "alpha_alpha" "nu_nu"       "mu_sigma"   
#>  [6] "mu_alpha"    "mu_nu"       "sigma_alpha" "sigma_nu"    "alpha_nu"   

# Against numerical differentiation of the log-likelihood: one component
# of the closed-form block, and two involving nu.
f <- function(v) sum(distrib_pdf(d, y, list(mu = v[1], sigma = v[2],
                                            alpha = v[3], nu = v[4]),
                                 log = TRUE))
H <- numDeriv::hessian(f, c(0, 1, 3, 6))
rbind(analytic = c(sum(h$mu_mu), sum(h$mu_nu), sum(h$nu_nu)),
      numeric = c(H[1, 1], H[1, 4], H[4, 4]))
#>               [,1]       [,2]       [,3]
#> analytic -6.590151 -0.4006573 0.04158258
#> numeric  -6.590170 -0.4006573 0.04158258

# The curvature in the location turns positive far out, as a Student t's
# does and a Gaussian's does not.
range(distrib_hess_y(d, seq(-20, 20, by = 0.5), th))
#> [1] -7.391840  1.643761
```
