# Multivariate Student t Observed Hessian

Computes the second derivatives of the log-density in closed form, by
differentiating the score of
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
once more. Every block picks up a term in \\\partial c/\partial\cdot\\,
the weight \\c = (\nu+p)/(\nu+q)\\ depending on the observation through
\\q\\. That dependence is what separates the family from the gaussian,
where \\c\\ is one and those terms are absent. Writing \\w =
\Sigma^{-1}(y-\mu)\\ and \\d = \nu + q\\, the location block is
\$\$\ell^{(\mu_a\mu_b)} = \frac{2c\\w_a w_b}{d} -
c\\(\Sigma^{-1})\_{ab},\$\$ and \\\partial c/\partial\nu = (q-p)/d^2\\
is what carries \\\nu\\ into every mixed block.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two differ here, `nu` carrying a log link by
  default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, one per unordered pair
of parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Unlike the gaussian's, no block here is constant across rows.

## Details

The pure \\\nu\\ component has three leading terms that cancel as
\\\nu\\ grows. They collapse exactly into
[`mvt_T()`](https://statmodels7.github.io/distributions7/reference/mvt_T.md),
with no series and no crossover for an even dimension; see
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
for what the direct form costs.

The quantity \\w^\top A_k w\\ appears in three of the blocks and is
formed once per free value, as is \\\Sigma^{-1}A_k\\.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\\eta\\ the free vector of the matrix
parametrization, \\A_k\\ and \\A\_{kl}\\ its first and second derivative
arrays, \\q\\ the squared Mahalanobis distance, \\c\\ the weight and
\\\ell^{(ij)}\\ the second derivative of the log-density in parameters
\\i\\ and \\j\\.

## See also

[`distrib_expected_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md)
for the expectation of this matrix, which is closed form,
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
for the score,
[`mvt_T()`](https://statmodels7.github.io/distributions7/reference/mvt_T.md)
for the cancellation in the \\\nu\\ component, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 25, theta)

H <- distrib_hessian(d, y, theta)
vapply(H, sum, numeric(1))
#>                   mu1_mu1                   mu2_mu2 sigma_log_L1_sigma_log_L1 
#>              -22.81097488              -35.16164546              -37.50229398 
#> sigma_log_L2_sigma_log_L2     sigma_L2.1_sigma_L2.1                     nu_nu 
#>              -23.57826977              -29.23625563               -0.08531611 
#>                   mu1_mu2          mu1_sigma_log_L1          mu1_sigma_log_L2 
#>               13.21985046              -11.66479639                0.71872483 
#>            mu1_sigma_L2.1                    mu1_nu          mu2_sigma_log_L1 
#>                3.68817570               -0.05242484                4.59023454 
#>          mu2_sigma_log_L2            mu2_sigma_L2.1                    mu2_nu 
#>               -2.07866833              -10.08058330               -0.06245101 
#> sigma_log_L1_sigma_log_L2   sigma_log_L1_sigma_L2.1           sigma_log_L1_nu 
#>                4.57647430                8.82634181                0.32240815 
#>   sigma_log_L2_sigma_L2.1           sigma_log_L2_nu             sigma_L2.1_nu 
#>               -2.70003054                0.25235942                0.01410377 

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
#> [1] 1.978684e-09

# The location block moves with the observation, where a gaussian's does
# not: the curvature at a far row is smaller.
z <- distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, theta))
near <- which.min(z$q); far <- which.max(z$q)
c(q_near = z$q[near], curv_near = H$mu1_mu1[near],
  q_far = z$q[far], curv_far = H$mu1_mu1[far])
#>      q_near   curv_near       q_far    curv_far 
#>  0.05803935 -1.33835717  9.77285228  0.10111337 
```
