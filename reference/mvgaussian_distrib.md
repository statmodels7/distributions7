# Construct a Multivariate Gaussian Distribution

The gaussian distribution on \\\mathbb{R}^p\\, with the mean a vector of
\\p\\ free parameters and the matrix carried by a structure from
parameters7.

## Usage

``` r
mvgaussian_distrib(n_dim, sigma = NULL, omega = NULL)
```

## Arguments

- n_dim:

  The dimension \\p\\.

- sigma:

  A parameters7 structure for the covariance. Defaults to
  `parameters7::log_cholesky(n_dim)` when neither structure is given.

- omega:

  A parameters7 structure for the precision.

## Value

An object of class
[`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md).

## Details

Exactly one of `sigma` and `omega` may be given, and the name of the
argument decides which side of the model the matrix parameter
parametrises: the covariance in the first case, the precision in the
second. One constructor returns one of two behaviours, in the manner of
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md),
which chooses between its continuous and discrete classes from the
arguments it is handed.

The precision form is the cheaper one and is worth preferring where the
modelling allows it. Written in \\\Omega\\, the log-density, the score
and the Hessian are multiplications, and the first term of the score is
the parameter's own `param_dlogdet()`; written in \\\Sigma\\ the same
quantities need a solve at every step.

**Parameters.** The mean contributes `mu1`, ..., `mup`, and the matrix
parameter contributes its free values prefixed by the matrix they
describe: `sigma_` for a covariance and `omega_` for a precision. A
two-dimensional gaussian on an unstructured covariance therefore has
five parameters, `mu1`, `mu2`, `sigma_log_L1`, `sigma_log_L2` and
`sigma_L2.1`, while the same structure on the precision gives
`omega_log_L1` and the rest. The prefix is what distinguishes the two
models in a printed table, since the name of a free value says how the
matrix is built and not which matrix it is.

All of the parameters are unconstrained, and their links are therefore
the identity: the constraint that makes the matrix positive definite
lives inside the matrix parameter, which is why it needs no link to
express it. A consequence worth knowing is that the parameter scale and
the link scale coincide here, so `scale = "link"` changes nothing.

**Reading a fit.** The free values are coordinates, not quantities
anybody reads.
[`mv_summary`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
carries the fit's variance matrix onto the standard deviations and
correlations by the delta method, and
[`print()`](https://rdrr.io/r/base/print.html) shows them; a precision
parametrisation also reports the conditional variances and the partial
correlations, which are what it describes directly. The conditional
variance is \\1/\Omega\_{jj} = \mathrm{Var}(Y_j \mid Y\_{-j})\\, and its
ratio to the marginal variance is \\1 - R_j^2\\ for the regression of
that coordinate on all the others.

**Rank.** A rank-deficient structure is refused. A singular covariance
gives a law supported on a subspace, with no density against Lebesgue
measure, and a singular precision gives a quadratic form that is flat
along its null space and does not normalise. The two are different
failures and both are failures; a structure of that kind is a legitimate
penalty and not a legitimate density.

**The response** is an \\n \times p\\ matrix, one row per observation. A
plain vector of length \\p\\ is read as a single observation.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)

## Examples

``` r
d <- mvgaussian_distrib(2)
d
#> Distribution: Multivariate Gaussian [2d, Sigma=log_cholesky]
#> Type:         Continuous, 2-dimensional
#> Dimensions:   multivariate
#> 
#> Parameters:
#>   mu1          (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   mu2          (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_log_L1 (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_log_L2 (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_L2.1   (covariance)         | Link: identity   | Domain: (-Inf, Inf)

theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)
y <- rbind(c(0, 0), c(1, -1))
distrib_pdf(d, y, theta, log = TRUE)
#> [1] -1.837877 -3.462877

# the covariance the free values describe
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25

# a diagonal covariance: two variances instead of three free values
mvgaussian_distrib(2, sigma = parameters7::diagonal_matrix(2))@params
#> [1] "mu1"          "mu2"          "sigma_log_d1" "sigma_log_d2"

# or the precision, which is the cheaper parametrisation
mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
```
