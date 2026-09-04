# Multivariate Gaussian Score

Computes the first derivatives of the log-density in closed form. With
\\w = \Sigma^{-1}(y - \mu)\\ and \\A_k\\ the derivative of \\\Sigma\\ in
the \\k\\-th free value of the matrix parametrization,
\$\$\frac{\partial \ell}{\partial \mu} = w, \qquad \frac{\partial
\ell}{\partial \eta_k} = -\frac{1}{2}\frac{\partial
\log\lvert\Sigma\rvert}{\partial \eta_k} + \frac{1}{2}\\ w^\top A_k
w.\$\$ The first term of the matrix component is the parametrization's
own `param_dlogdet()`, so no trace is formed here; its sign is flipped
where the parametrization carries the precision. The mean component
needs no derivative array at all, which is why the mean and the matrix
cost so differently in a fit.

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
  scales coincide and the results are the same numbers to the bit.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector of length \\n\\ per parameter, in
`distrib@params` order: \\p\\ mean components followed by one per free
value of the matrix parametrization.

## Notation

\\\mu\\ is the mean, \\\Sigma\\ the covariance, \\\eta\\ the free vector
of the matrix parametrization, \\A_k = \partial\Sigma/\partial\eta_k\\
and \\\ell\\ the log-density of one observation.

## See also

[`distrib_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
for the observed curvature,
[`distrib_expected_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md)
for the information,
[`distrib_grad_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvGaussianDistrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 25, theta)

g <- distrib_gradient(d, y, theta)
names(g)
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
vapply(g, sum, numeric(1))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1 
#>    3.4591533    0.9841865   -2.1249520  -13.0020520   -1.2473141 

# Against a numerical derivative of the log-likelihood, which shares no
# arithmetic with the closed form.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 4.368311e-09

# The mean component is the whitened residual, one row per observation.
S <- mv_sigma(d, theta)
w <- sweep(y, 2, c(0.5, -0.3)) %*% solve(S)
all.equal(cbind(g$mu1, g$mu2), w, check.attributes = FALSE)
#> [1] TRUE

# Every link is the identity, so the link scale changes nothing.
identical(g, distrib_gradient(d, y, theta, scale = "link"))
#> [1] TRUE
```
