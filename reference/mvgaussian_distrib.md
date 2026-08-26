# Construct a Multivariate Gaussian Distribution

Builds the gaussian family on \\\mathbb{R}^p\\. The mean is a vector of
\\p\\ free parameters and the matrix is carried by a parameters7
parametrization, either as the covariance \\\Sigma\\ or as the precision
\\\Omega = \Sigma^{-1}\\. The free values of that parametrization become
ordinary scalar parameters of the distribution, so the object answers
every generic of the `distrib` contract with `theta` a flat named list
of numbers. The default is an unstructured covariance in the
log-Cholesky parametrization, which is `p * (p + 1) / 2` free values.

## Usage

``` r
mvgaussian_distrib(n_dim, sigma = NULL, omega = NULL)
```

## Arguments

- n_dim:

  The dimension \\p\\ of one observation. A single positive whole
  number, finite and at least 1. Anything else throws an error.

- sigma:

  A parameters7 parametrization of the covariance, of dimension `n_dim`
  and of full rank. Defaults to `parameters7::log_cholesky(n_dim)` when
  neither this nor `omega` is given.

- omega:

  A parameters7 parametrization of the precision, of dimension `n_dim`
  and of full rank. Defaults to `NULL`. Giving both this and `sigma` is
  an error.

## Value

An S7 object of class
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md),
with `param` the parametrization supplied and `inverted` recording which
side it carries. Its `params` are `mu1`, ..., `mup` followed by the
prefixed free names, `n_params` is `p + param@n_free`, every entry of
`params_bounds` is \\(-\infty, \infty)\\ and every link is the identity.

## Which side the matrix parametrizes

Give at most one of `sigma` and `omega`, and the name of the argument
decides which side of the model the matrix parametrization describes.
Giving both is an error: the two name the same matrix and a distribution
parametrized by both would be over-determined.

The two forms describe the same family and cost about the same. Every
derivative here is written in the covariance, so a precision
parametrization is inverted once per call and its derivative arrays are
carried across by \\\partial\Sigma/\partial\eta_k = -\Sigma A_k
\Sigma\\. Measured at \\n = 1000\\ on an unstructured matrix, the
precision form costs 0.86 to 1.07 times the covariance form for the
score and 0.97 to 1.10 for the Hessian, across \\p = 4, 8, 12\\. Choose
the side the model is written in.

## Parameters

The mean contributes `mu1`, ..., `mup`. The matrix contributes the
parametrization's free values, prefixed by the matrix they describe:
`sigma_` for a covariance and `omega_` for a precision. A
two-dimensional gaussian on an unstructured covariance therefore has the
five parameters `mu1`, `mu2`, `sigma_log_L1`, `sigma_log_L2` and
`sigma_L2.1`, while the same parametrization on the precision gives
`omega_log_L1` and the rest. The prefix is what distinguishes the two
models in a printed table: a free value's name says how the matrix is
built, not which matrix it is.

All of these parameters are unconstrained, so all of the links are the
identity. The constraint that makes the matrix positive definite lives
inside the matrix parametrization, which needs no link to express it.
One practical consequence is that the parameter scale and the link scale
coincide, so `scale = "link"` returns the same numbers as
`scale = "parameter"`, to the bit.

## Reading a fit

The free values are coordinates, and nobody reads a logarithm of a
diagonal entry of a Cholesky factor.
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
carries a fit's variance matrix onto the standard deviations and the
correlations by the delta method, and
[`print()`](https://rdrr.io/r/base/print.html) shows those. A precision
parametrization also reports the conditional standard deviations and the
partial correlations, which are what it describes directly:
\\1/\Omega\_{jj} = \operatorname{Var}(Y_j \mid Y\_{-j})\\, whose ratio
to the marginal variance is \\1 - R_j^2\\ for the regression of that
coordinate on all the others. At \\p = 2\\ the partial correlation and
the marginal correlation are the same number, so it is not printed
twice.

## Rank

A rank-deficient parametrization is rejected with an error. A singular
covariance gives a law supported on a subspace, with no density against
Lebesgue measure; a singular precision gives a quadratic form that is
flat along its null space and does not normalize. Both are failures, and
a parametrization of that kind is a legitimate penalty but not a
legitimate density.

## The response

An \\n \times p\\ matrix, one row per observation. A plain vector of
length \\p\\ is read as a single observation. A parameter may not vary
from observation to observation here: the parametrization describes one
matrix for the whole sample.

## The distribution

\$\$f(y) = (2\pi)^{-p/2}\lvert \Sigma
\rvert^{-1/2}\exp\\\left\\-\tfrac{1}{2}(y-\mu)^\top\Sigma^{-1}(y-\mu)\right\\\$\$
on \\y \in \mathbb{R}^{p}\\, with \$\$\mathbb{E}\[Y\] = \mu, \qquad
\operatorname{Var}(Y) = \Sigma.\$\$

## Notation

\\\mu\\ is the mean vector, \\\Sigma\\ the covariance, \\\Omega\\ the
precision, \\p\\ the dimension of one observation and \\n\\ the number
of observations. \\\eta\\ is the free vector of the matrix
parametrization, the unconstrained coordinates an optimizer moves, and
\\A_k = \partial\Sigma/\partial\eta_k\\.

## See also

[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the one-dimensional case,
[`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
for the heavy-tailed sibling,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the quantities a fit reports,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate one, and
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
for the default parametrization.

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

theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)

# The covariance the free values describe, and its positive definiteness,
# which holds for any five numbers at all.
S <- mv_sigma(d, theta)
round(S, 4)
#>        v1     v2
#> v1 1.2214 0.4421
#> v2 0.4421 0.8303
eigen(S, only.values = TRUE)$values
#> [1] 1.5092462 0.5424766

# The density, against the formula written out.
y <- rbind(c(0, 0), c(1, -1))
mu <- c(0.5, -0.3)
ref <- apply(y, 1, function(r)
  -0.5 * (2 * log(2 * pi) + log(det(S)) + t(r - mu) %*% solve(S) %*% (r - mu)))
all.equal(distrib_pdf(d, y, theta, log = TRUE), as.numeric(ref))
#> [1] TRUE

# A structured covariance spends fewer parameters: two variances, no
# correlation.
mvgaussian_distrib(2, sigma = parameters7::diagonal_matrix(2))@params
#> [1] "mu1"          "mu2"          "sigma_log_d1" "sigma_log_d2"

# The precision form carries the same law and reports its own quantities.
o <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))
o@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
             omega_log_L2 = -0.2, omega_L2.1 = 0.4)
Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
all.equal(mv_sigma(o, th_o), solve(Om), check.attributes = FALSE)
#> [1] TRUE

# The conditional variance of the first coordinate given the second, and
# the marginal variance beside it.
c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#> conditional    marginal 
#>   0.8187308   1.0141552 
```
