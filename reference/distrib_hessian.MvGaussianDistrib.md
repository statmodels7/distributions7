# Multivariate Gaussian Observed Hessian

Computes the second derivatives of the log-density in closed form. With
\\w = \Sigma^{-1}(y-\mu)\\ and \\A_k\\, \\A\_{kl}\\ the first and second
derivatives of \\\Sigma\\ in the free values, \$\$\ell^{(\mu_a \mu_b)} =
-(\Sigma^{-1})\_{ab}, \qquad \ell^{(\mu_a \eta_k)} = -(\Sigma^{-1} A_k
w)\_a,\$\$ \$\$\ell^{(\eta_k \eta_l)} = -\tfrac{1}{2}\frac{\partial^2
\log\lvert\Sigma\rvert}{\partial\eta_k \partial\eta_l} + \tfrac{1}{2}\\
w^\top A\_{kl} w - w^\top A_l \Sigma^{-1} A_k w.\$\$ The mean block does
not depend on the response, so it is constant across rows; the mixed and
matrix blocks do. The product \\\Sigma^{-1}A_k\\ is formed once and
reused between the mixed block and the matrix block.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. Every link of this family is the identity, so the two
  scales coincide.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, one per unordered pair
of parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md):
the diagonal first, then the off-diagonal pairs.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\\eta\\ the free vector
of the matrix parametrization, \\A_k\\ and \\A\_{kl}\\ its first and
second derivative arrays, and \\\ell^{(ij)}\\ the second derivative of
the log-density in parameters \\i\\ and \\j\\.

## See also

[`distrib_expected_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md)
for the expectation of this matrix, which is simpler,
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
for the score, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 25, theta)

H <- distrib_hessian(d, y, theta)
names(H)
#>  [1] "mu1_mu1"                   "mu2_mu2"                  
#>  [3] "sigma_log_L1_sigma_log_L1" "sigma_log_L2_sigma_log_L2"
#>  [5] "sigma_L2.1_sigma_L2.1"     "mu1_mu2"                  
#>  [7] "mu1_sigma_log_L1"          "mu1_sigma_log_L2"         
#>  [9] "mu1_sigma_L2.1"            "mu2_sigma_log_L1"         
#> [11] "mu2_sigma_log_L2"          "mu2_sigma_L2.1"           
#> [13] "sigma_log_L1_sigma_log_L2" "sigma_log_L1_sigma_L2.1"  
#> [15] "sigma_log_L2_sigma_L2.1"  

# The mean block is minus the inverse covariance and does not move with y.
c(H$mu1_mu1[1], H$mu1_mu2[1])
#> [1] -1.0141552  0.5399435
-solve(mv_sigma(d, theta))[1, ]
#>         v1         v2 
#> -1.0141552  0.5399435 

# Against a numerical Hessian of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
Hn <- numDeriv::hessian(ll, unlist(theta))
ref <- vapply(distributions7:::hess_pairs(d@params),
              function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
              numeric(1))
max(abs(vapply(H, sum, numeric(1)) - ref))
#> [1] 1.013385e-08
```
