# Construct a Multivariate Gaussian Distribution

The gaussian distribution on \\\mathbb{R}^p\\, with the mean a vector of
\\p\\ free parameters and the matrix carried by a structure from
covstructs7.

## Usage

``` r
mvgaussian_distrib(n_dim, struct_sigma = NULL, struct_omega = NULL)
```

## Arguments

- n_dim:

  The dimension \\p\\.

- struct_sigma:

  A covstructs7 structure for the covariance. Defaults to
  `covstructs7::log_cholesky(n_dim)` when neither structure is given.

- struct_omega:

  A covstructs7 structure for the precision.

## Value

An object of class
[`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md).

## Details

Exactly one of `struct_sigma` and `struct_omega` may be given, and the
name of the argument decides which side of the model the structure
parametrises: the covariance in the first case, the precision in the
second. One constructor returns one of two behaviours, in the manner of
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md),
which chooses between its continuous and discrete classes from the
arguments it is handed.

The precision form is the cheaper one and is worth preferring where the
modelling allows it. Written in \\\Omega\\, the log-density, the score
and the Hessian are multiplications, and the first term of the score is
the structure's own `struct_dlogdet()`; written in \\\Sigma\\ the same
quantities need a solve at every step.

**Parameters.** The mean contributes `mu1`, ..., `mup` and the structure
contributes its own free values under their own names, so a
two-dimensional gaussian on an unstructured covariance has five
parameters: `mu1`, `mu2`, `log_L1`, `log_L2`, `L2.1`. All of them are
unconstrained, and their links are therefore the identity: the
constraint that makes the matrix positive definite lives inside the
structure, which is why it needs no link to express it. A consequence
worth knowing is that the parameter scale and the link scale coincide
here, so `scale = "link"` changes nothing.

**Rank.** A rank-deficient structure is refused. A singular covariance
gives a law supported on a subspace, with no density against Lebesgue
measure, and a singular precision gives a quadratic form that is flat
along its null space and does not normalise. The two are different
failures and both are failures; a structure of that kind is a legitimate
penalty and not a legitimate density.

**The response** is an \\n \times p\\ matrix, one row per observation. A
plain vector of length \\p\\ is read as a single observation.

## See also

[`gaussian_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md),
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`log_cholesky`](https://statmodels7.github.io/covstructs7/reference/log_cholesky.html)

## Examples

``` r
d <- mvgaussian_distrib(2)
d
#> Distribution: Multivariate Gaussian [2d, Sigma=log_cholesky]
#> Type:         Continuous, 2-dimensional
#> Dimensions:   multivariate
#> 
#> Parameters:
#>   mu1    (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   mu2    (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   log_L1 (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#>   log_L2 (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#>   L2.1   (covariance)         | Link: identity   | Domain: (-Inf, Inf)

theta <- list(mu1 = 0, mu2 = 0, log_L1 = 0, log_L2 = 0, L2.1 = 0.5)
y <- rbind(c(0, 0), c(1, -1))
distrib_pdf(d, y, theta, log = TRUE)
#> [1] -1.837877 -3.462877

# the covariance the free values describe
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25

# a diagonal covariance: two variances instead of three free values
mvgaussian_distrib(2, struct_sigma = covstructs7::diag_struct(2))@params
#> [1] "mu1" "mu2" "d1"  "d2" 

# or the precision, which is the cheaper parametrisation
mvgaussian_distrib(2, struct_omega = covstructs7::log_cholesky(2))@params
#> [1] "mu1"    "mu2"    "log_L1" "log_L2" "L2.1"  
```
