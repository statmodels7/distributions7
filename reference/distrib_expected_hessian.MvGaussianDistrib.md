# Multivariate Gaussian Expected Information

Computes the expectation of the observed Hessian in closed form, which
for this family is simpler than the observed matrix itself:
\$\$\mathbb{E}\[\ell^{(\mu_a \mu_b)}\] = -(\Sigma^{-1})\_{ab}, \qquad
\mathbb{E}\[\ell^{(\mu_a \eta_k)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\eta_k \eta_l)}\] =
-\tfrac{1}{2}\operatorname{tr}(\Sigma^{-1} A_k \Sigma^{-1} A_l).\$\$ No
component depends on the data, so every returned vector is constant
across rows and `y` is read for its row count alone.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. Only its row count
  is used: the expectation is taken over the law, so no observation
  enters any component.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. Every link of this family is the identity, so the two
  scales coincide.

- approx:

  Ignored: the expectation is exact and no approximation strategy is
  consulted. Present so that the signature matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each vector constant.

## Why the mixed block vanishes

Each mixed component is a linear function of \\w = \Sigma^{-1}(y-\mu)\\,
and \\\mathbb{E}\[w\] = 0\\. The mean parameters and the matrix
parameters are therefore orthogonal, and the information matrix is block
diagonal whatever the matrix parametrization is. Fisher scoring on this
family is well behaved for that reason: a step in the mean and a step in
the covariance do not interfere.

## Why no second derivative array is needed

The observed matrix block carries \\A\_{kl}\\ through both the
log-determinant term and the quadratic form. Under the expectation the
two cancel, and only the first derivatives survive. That saves the whole
`param_d2()` computation, which is the dearest part of the observed
Hessian at any dimension worth the name.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\\eta\\ the free vector
of the matrix parametrization, \\A_k = \partial\Sigma/\partial\eta_k\\,
and \\\ell^{(ij)}\\ the second derivative of the log-density in
parameters \\i\\ and \\j\\.

## See also

[`distrib_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
for the observed matrix,
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
for the score,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
whose Fisher scoring inverts this matrix, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
y <- matrix(0, 3, 2)

EH <- distrib_expected_hessian(d, y, theta)

# The mean block is minus the inverse covariance.
c(EH$mu1_mu1[1], EH$mu1_mu2[1])
#> [1] -1.0141552  0.5399435
-solve(mv_sigma(d, theta))[1, ]
#>         v1         v2 
#> -1.0141552  0.5399435 

# Every mixed mean-matrix component is exactly zero, not merely small.
mixed <- grep("^mu[0-9]+_sigma", names(EH), value = TRUE)
vapply(EH[mixed], function(z) z[1], numeric(1))
#> mu1_sigma_log_L1 mu1_sigma_log_L2   mu1_sigma_L2.1 mu2_sigma_log_L1 
#>                0                0                0                0 
#> mu2_sigma_log_L2   mu2_sigma_L2.1 
#>                0                0 

# And the closed form is what averaging the observed Hessian converges to.
set.seed(4)
big <- distrib_rng(d, 50000, theta)
emp <- vapply(distrib_hessian(d, big, theta), mean, numeric(1))
round(rbind(sampled = emp,
            closed = vapply(EH, function(z) z[1], numeric(1))), 3)
#>         mu1_mu1 mu2_mu2 sigma_log_L1_sigma_log_L1 sigma_log_L2_sigma_log_L2
#> sampled  -1.014  -1.492                    -2.224                        -2
#> closed   -1.014  -1.492                    -2.239                        -2
#>         sigma_L2.1_sigma_L2.1 mu1_mu2 mu1_sigma_log_L1 mu1_sigma_log_L2
#> sampled                -1.482    0.54           -0.006           -0.003
#> closed                 -1.492    0.54            0.000            0.000
#>         mu1_sigma_L2.1 mu2_sigma_log_L1 mu2_sigma_log_L2 mu2_sigma_L2.1
#> sampled          0.005            0.002            0.007         -0.004
#> closed           0.000            0.000            0.000          0.000
#>         sigma_log_L1_sigma_log_L2 sigma_log_L1_sigma_L2.1
#> sampled                    -0.001                   0.593
#> closed                      0.000                   0.597
#>         sigma_log_L2_sigma_L2.1
#> sampled                   0.001
#> closed                    0.000
```
