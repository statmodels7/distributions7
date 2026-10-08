# NB1 Expected Hessian

Returns the expectation of the observed Hessian under the model. Every
term carrying \\P = \psi(y+r) - \psi(r) - \log(1+\theta)\\ drops out,
its expectation vanishing by the first Bartlett identity. With \\u_j =
\mu + \theta j\\, the entries are \\-A\\, \\W/\theta^2\\ and
\\-W/(\mu\theta)\\, where \$\$A = \mathbb{E}\Big\[\sum\_{j\<Y}
u_j^{-2}\Big\], \qquad W = \mathbb{E}\Big\[\sum\_{j\<Y}
\Big(\frac{1}{1+\theta} - \frac{\mu^2}{u_j^2}\Big)\Big\].\$\$ Neither
has a closed form. Each is summed against the exact mass until the terms
fall below a tolerance, so this is a truncated exact sum and not a
quadrature or a simulation. `approx` and `nsim` are ignored, and `y` is
read only for its length.

**The mixed entry does not vanish**, so the mean and the dispersion are
not orthogonal in this family. Measured at four settings it is 0.0172,
0.0364, 0.0108 and 0.0157; in
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
the same entry is exactly zero. That is a difference between the two
negative binomials rather than a difference of parametrization.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- y:

  A numeric vector of counts. Only its length is used.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored. The surviving expectation is a truncated exact sum. Accepted
  so that the signature matches the generic's.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_theta` and
`theta_theta`, in that order, each of length
`max(length(y), length(mu), length(theta))` and constant within itself
when the parameters are.

## Accuracy at small theta

The summand of \\W\\ is a rational function of \\(\mu, \theta)\\,
\\\theta(2\mu j + \theta j^2 - \mu^2)/((1+\theta) u_j^2)\\, so the terms
of order \\\mu/\theta^2\\ that a composition through the size \\r =
\mu/\theta\\ forms are never formed. The dispersion entry still loses
digits as \\\theta \to 0\\, where the summands change sign and their sum
tends to a finite limit. Measured at \\\mu = 4\\, it reads
\\-0.49999887\\ at \\\theta = 10^{-6}\\, \\-0.49999994\\ at \\10^{-8}\\
and \\-0.5000013\\ at \\10^{-10}\\, against the limit \\-1/2\\, and the
matrix stays negative definite.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
first Bartlett identity is \\\mathbb{E}\[\partial\ell/\partial\theta\] =
0\\.

## See also

[`distrib_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md)
for the quadratic-variance family, whose mixed entry is exactly zero,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
th <- list(mu = 4, theta = 4)
e <- distrib_expected_hessian(d, c(0, 2, 6), th)
lapply(e, unique)
#> $mu_mu
#> [1] -0.06717466
#> 
#> $mu_theta
#> [1] 0.01717466
#> 
#> $theta_theta
#> [1] -0.01717466
#> 

# Negative definite at a moderate dispersion.
M <- matrix(c(e$mu_mu[1], e$mu_theta[1], e$mu_theta[1], e$theta_theta[1]), 2)
eigen(M, only.values = TRUE)$values
#> [1] -0.01184367 -0.07250565

# The mixed entry is not zero, so the mean and the dispersion are not
# orthogonal here; in the quadratic-variance family it is exactly zero.
c(nb1 = e$mu_theta[1],
  nb2 = distrib_expected_hessian(negbin2_distrib(), 0, th)$mu_theta)
#>        nb1        nb2 
#> 0.01717466 0.00000000 

# The dispersion entry loses its digits as theta goes to zero, and turns
# positive, which an expected second derivative cannot be.
vapply(c(1e-2, 1e-4, 1e-6, 1e-8),
       function(t) distrib_expected_hessian(d, 0,
                     list(mu = 4, theta = t))$theta_theta,
       numeric(1))
#> [1] -0.4889418 -0.4998875 -0.4999989 -0.4999999
```
