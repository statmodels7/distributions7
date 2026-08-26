# Student's t Mixed Derivatives

Closed form in all three parameters. With \\r = y - \mu\\ and \\D =
\nu\sigma^2 + r^2\\, \$\$\frac{\partial^2 \ell}{\partial y\\ \partial
\mu} = \frac{(\nu+1)(\nu\sigma^2 - r^2)}{D^2}, \qquad \frac{\partial^2
\ell}{\partial y\\ \partial \sigma} = \frac{2\nu\sigma(\nu+1) r}{D^2},
\qquad \frac{\partial^2 \ell}{\partial y\\ \partial \nu} =
-\frac{r(r^2 - \sigma^2)}{D^2}.\$\$

The \\\nu\\ component is elementary here, unlike the \\\nu\\ components
of the cdf derivatives: the log-density carries `lgamma` and a logarithm
of \\D\\, both differentiable in \\\nu\\ in closed form, and no
distribution function appears.

Note that the \\\mu\\ component changes sign at \\\|r\| =
\sigma\sqrt\nu\\ and decays as \\r^{-2}\\, where the Gaussian's is
constant. That is the redescending score seen one derivative along.

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, aligned by the
  generic. Each may be length 1 or `length(y)`.

- scale:

  Handled by the generic after dispatch; this method always returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list with components `mu`, `sigma` and `nu`, each a numeric
vector of length `length(y)`.

## See also

[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic;
[`distrib_cross_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gaussian1Distrib.md),
the limit as \\\nu\\ grows;
[`distrib_grad_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.StudentT1Distrib.md)
for the redescending score itself.

## Examples

``` r
d <- student_t1_distrib()
th <- list(mu = 0.2, sigma = 1.1, nu = 6)
distrib_cross_y(d, c(-1, 0, 2), th)
#> $mu
#> [1] 0.5382481 0.9483956 0.2552381
#> 
#> $sigma
#> [1] -1.4649227 -0.3467818  1.5085714
#> 
#> $nu
#> [1]  0.003646453 -0.004391068 -0.033142857
#> 

# The mu component changes sign where the residual passes sigma*sqrt(nu),
# which is where the score of a t stops growing and starts to redescend.
r0 <- 1.1 * sqrt(6)
distrib_cross_y(d, 0.2 + c(0.5, 1, 1.5) * r0, th)$mu
#> [1]  4.628099e-01  5.897868e-17 -1.141050e-01
```
