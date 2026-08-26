# Student's t Mixed Second-Response Derivatives

Closed form in all three parameters. With \\z = (y-\mu)/\sigma\\ and
\\Q(z) = (\nu - z^2)/(\nu + z^2)^2\\ the response curvature is
\\\ell^{(yy)} = -(\nu+1)Q/\sigma^2\\, so writing \\Q' = 2z(z^2 -
3\nu)/(\nu+z^2)^3\\ and \\\partial Q/\partial\nu =
(3z^2-\nu)/(\nu+z^2)^3\\, \$\$\frac{\partial\ell^{(yy)}}{\partial\mu} =
\frac{(\nu+1)Q'}{\sigma^3}, \qquad
\frac{\partial\ell^{(yy)}}{\partial\sigma} = \frac{(\nu+1)(2Q +
zQ')}{\sigma^3}, \qquad \frac{\partial\ell^{(yy)}}{\partial\nu} =
-\frac{Q + (\nu+1)\partial\_\nu Q}{\sigma^2}.\$\$

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with `mu`, `sigma` and `nu`.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors of length `length(y)`, keyed `mu`,
`sigma` and `nu`.

## Details

Nothing vanishes here, unlike the gaussian's, and every component varies
with the data: the t's response curvature is a function of \\z\\ rather
than a constant, and it is that dependence which makes the family
redescending. The three formulas share \\Q\\, \\Q'\\ and \\\partial\_\nu
Q\\, so each is evaluated once per observation.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other. \\\ell^{(yy)}\\ is \\\partial^2\ell/\partial y^2\\.

## See also

[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic,
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the quantity being differentiated, and
[`distrib_cross2_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.Gaussian1Distrib.md),
the limit as \\\nu\\ grows.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3, nu = 6)
distrib_cross2_y(d, y, theta)
#> $mu
#> [1]  0.3148726  0.1550762 -0.3046618
#> 
#> $sigma
#> [1] 0.2625888 0.9653573 0.1311434
#> 
#> $nu
#> [1] -0.02742339  0.01050930 -0.03278889
#> 

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2],
                                              nu = v[3]))
numDeriv::grad(f, c(0.4, 1.3, 6))
#> [1]  0.31487259  0.26258876 -0.02742339

# Every component varies with y, where the gaussian's do not.
distrib_cross2_y(d, y, theta)$mu
#> [1]  0.3148726  0.1550762 -0.3046618
distrib_cross2_y(gaussian1_distrib(), y, theta[1:2])$mu
#> [1] 0 0 0
```
