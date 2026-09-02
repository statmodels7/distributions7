# Multivariate Student t Expected Information

Computes the expectation of the observed Hessian in closed form, taken
from the family's own scale mixture. Every component is a Beta moment or
a polygamma, so no sampling and no quadrature runs, and two calls at the
same parameters return the same numbers. No component depends on the
data, so every returned vector is constant across rows and `y` is read
for its row count alone.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. Only its row count
  is used: the expectation is taken over the law, so no observation
  enters any component.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two differ here, `nu` carrying a log link by
  default.

- approx:

  Ignored: the expectation is exact and no approximation strategy is
  consulted.
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  rejects a `fisher_scoring(approx =)` for this family for that reason.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each vector constant. The returned quantities are the expected SECOND
DERIVATIVES, so the information is their negative.

## Why the expectations separate

Writing \\q = z^\top\Sigma^{-1}z\\ and \\c = (\nu+p)/(\nu+q)\\, the
variable \\v = q/(q+\nu)\\ is EXACTLY \\\mathrm{Beta}(p/2, \nu/2)\\ and
is independent of the direction \\z/\lVert z\rVert\\, which is uniform
on the sphere. Each expectation therefore factors into a radial part and
an angular part: the radial one gives \\cq = (\nu+p)v\\ and \\c^2q^2 =
(\nu+p)^2v^2\\, and the angular one gives \\\mathbb{E}\[ee^\top\] =
I/p\\ and \\\mathbb{E}\[(e^\top Be)(e^\top Ce)\] =
\\\operatorname{tr}B\operatorname{tr}C + 2\operatorname{tr}(BC)\\/
\\p(p+2)\\\\.

## The blocks

With \\t_k = \operatorname{tr}(\Sigma^{-1}A_k)\\ and \\T\_{kl} =
\operatorname{tr}(\Sigma^{-1}A_k\Sigma^{-1}A_l)\\, \$\$I\_{\mu\mu} =
\frac{\nu+p}{\nu+p+2}\\\Sigma^{-1}, \qquad I\_{kl} =
\frac{(\nu+p)T\_{kl} - t_kt_l}{2(\nu+p+2)}, \qquad I\_{k\nu} =
\frac{-t_k}{(\nu+p)(\nu+p+2)},\$\$ \$\$I\_{\nu\nu} =
\frac{\psi'(\nu/2) - \psi'\\(\nu+p)/2\\}{4} - \frac{p}{\nu(\nu+p)} +
\frac{p}{2\nu(\nu+p+2)}.\$\$ The location is orthogonal to everything
else, every cross-expectation with it being odd in \\z\\.

## The gaussian limit

Each block reduces to the gaussian's as \\\nu \to \infty\\:
\\I\_{\mu\mu} \to \Sigma^{-1}\\, \\I\_{kl} \to T\_{kl}/2\\, and
\\I\_{\nu\nu} \to 0\\, the degrees of freedom ceasing to be identified
in the limit. Measured at \\p = 2\\ on an unstructured matrix,
\\I\_{\nu\nu}\\ runs \\2.8\times10^{-3}\\, \\8.2\times10^{-6}\\,
\\4.9\times10^{-9}\\, \\1.3\times10^{-14}\\ at \\\nu = 6, 30, 200,
5000\\, which is why a standard error for \\\nu\\ is worth reading only
where the tail is genuinely heavy.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\\eta\\ the free vector of the matrix
parametrization, \\A_k = \partial\Sigma/\partial\eta_k\\, \\\psi'\\ the
trigamma function, \\e\\ a uniform direction on the unit sphere and
\\I\\ an expected information block.

## See also

[`distrib_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md)
for the observed matrix,
[`distrib_expected_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md)
for the limiting family,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
whose Fisher scoring inverts this matrix, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

EH <- distrib_expected_hessian(d, matrix(0, 3, 2), theta)
round(vapply(EH, function(z) z[1], numeric(1)), 5)
#>                   mu1_mu1                   mu2_mu2 sigma_log_L1_sigma_log_L1 
#>                  -0.81132                  -1.19346                  -1.59095 
#> sigma_log_L2_sigma_log_L2     sigma_L2.1_sigma_L2.1                     nu_nu 
#>                  -1.40000                  -1.19346                  -0.00278 
#>                   mu1_mu2          mu1_sigma_log_L1          mu1_sigma_log_L2 
#>                   0.43195                   0.00000                   0.00000 
#>            mu1_sigma_L2.1                    mu1_nu          mu2_sigma_log_L1 
#>                   0.00000                   0.00000                   0.00000 
#>          mu2_sigma_log_L2            mu2_sigma_L2.1                    mu2_nu 
#>                   0.00000                   0.00000                   0.00000 
#> sigma_log_L1_sigma_log_L2   sigma_log_L1_sigma_L2.1           sigma_log_L1_nu 
#>                   0.20000                   0.47738                   0.02500 
#>   sigma_log_L2_sigma_L2.1           sigma_log_L2_nu             sigma_L2.1_nu 
#>                   0.00000                   0.02500                   0.00000 

# Closed form, so two calls agree to the bit and nothing is sampled.
identical(EH, distrib_expected_hessian(d, matrix(0, 3, 2), theta))
#> [1] TRUE

# Every location-matrix and location-nu component is exactly zero.
orth <- grep("^mu[0-9]+_(sigma|nu)", names(EH), value = TRUE)
vapply(EH[orth], function(z) z[1], numeric(1))
#> mu1_sigma_log_L1 mu1_sigma_log_L2   mu1_sigma_L2.1           mu1_nu 
#>                0                0                0                0 
#> mu2_sigma_log_L1 mu2_sigma_log_L2   mu2_sigma_L2.1           mu2_nu 
#>                0                0                0                0 

# And the closed form is what averaging the observed Hessian converges to.
set.seed(3)
big <- distrib_rng(d, 50000, theta)
round(rbind(sampled = vapply(distrib_hessian(d, big, theta), mean, numeric(1)),
            closed = vapply(EH, function(z) z[1], numeric(1))), 3)
#>         mu1_mu1 mu2_mu2 sigma_log_L1_sigma_log_L1 sigma_log_L2_sigma_log_L2
#> sampled  -0.810  -1.192                    -1.597                    -1.398
#> closed   -0.811  -1.193                    -1.591                    -1.400
#>         sigma_L2.1_sigma_L2.1  nu_nu mu1_mu2 mu1_sigma_log_L1 mu1_sigma_log_L2
#> sampled                -1.206 -0.003   0.433                0                0
#> closed                 -1.193 -0.003   0.432                0                0
#>         mu1_sigma_L2.1 mu1_nu mu2_sigma_log_L1 mu2_sigma_log_L2 mu2_sigma_L2.1
#> sampled              0      0            0.002           -0.003         -0.002
#> closed               0      0            0.000            0.000          0.000
#>         mu2_nu sigma_log_L1_sigma_log_L2 sigma_log_L1_sigma_L2.1
#> sampled      0                     0.204                   0.478
#> closed       0                     0.200                   0.477
#>         sigma_log_L1_nu sigma_log_L2_sigma_L2.1 sigma_log_L2_nu sigma_L2.1_nu
#> sampled           0.026                  -0.009           0.025             0
#> closed            0.025                   0.000           0.025             0

# The nu block vanishes as the tail lightens.
vapply(c(6, 30, 200, 5000), function(nu) {
  t2 <- theta; t2$nu <- nu
  -distrib_expected_hessian(d, matrix(0, 1, 2), t2)$nu_nu[1]
}, numeric(1))
#> [1] 2.777778e-03 8.169935e-06 4.853427e-09 1.278470e-14
```
