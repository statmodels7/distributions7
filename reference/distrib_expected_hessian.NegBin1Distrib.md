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

The mixed entry does not vanish, so the mean and the dispersion are not
orthogonal in this family; in
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
the same entry is exactly zero. This is a difference between the two
families, not between two parametrizations of one.

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
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- approx:

  Ignored. The surviving expectation is a truncated exact sum. It is
  accepted so that the signature matches the generic's.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. For
  fewer observations than an internal threshold the kernel runs
  sequentially whatever the count. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_theta` and
`theta_theta`, in that order, each of length
`max(length(y), length(mu), length(theta))` and constant within itself
when the parameters are.

## Accuracy at small theta

The summand of \\W\\ is a rational function of \\(\mu, \theta)\\,
\\\theta(2\mu j + \theta j^2 - \mu^2)/((1+\theta) u_j^2)\\, so the terms
of order \\\mu/\theta^2\\ that a composition through the size \\r =
\mu/\theta\\ forms are never formed. The summands still change sign as
\\\theta \to 0\\, and \\W/\theta^2\\ loses relative accuracy in
proportion to \\1/\theta\\, so \\W\\ is used only for \\\theta \ge
0.02\\. For \\\theta \< 0.02\\ the dispersion entry is computed as the
expectation of the observed one in the form without the size
([`distrib_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md)),
\$\$\mathbb{E}\[\ell^{(\theta\theta)}\] = -\mathbb{E}\Big\[\sum\_{j\<Y}
\frac{j^2}{u_j^2}\Big\] - \mu q'(\theta), \qquad q(\theta) =
\frac{\theta - \log(1+\theta)}{\theta^2},\$\$ whose summands are all
positive, and the mixed entry is \\-\theta/\mu\\ times it. Both keep
their accuracy down to the Poisson limit, where the dispersion entry
tends to \\-1/2\\. The derivatives of the expected information
([`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md))
are sums over the same mass and take the same two routes.

## Notation

The expected information is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the observed information under the model. The first
Bartlett identity is \\\mathbb{E}\[\partial\ell/\partial\theta\] = 0\\.

## See also

[`distrib_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md)
for the observed quantity whose expectation this is,
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

# As theta goes to zero the dispersion entry tends to its limit, -1/2.
vapply(c(1e-2, 1e-4, 1e-6, 1e-8),
       function(t) distrib_expected_hessian(d, 0,
                     list(mu = 4, theta = t))$theta_theta,
       numeric(1))
#> [1] -0.4889418 -0.4998875 -0.4999989 -0.5000000
```
