# Construct a Multivariate Gaussian Distribution

Two families on \\\mathbb{R}^p\\, differing in which matrix of the law
their parametrization describes: `mvgaussian1_distrib()` parametrizes
the covariance \\\Sigma\\ and `mvgaussian2_distrib()` the precision
\\\Omega = \Sigma^{-1}\\. The mean is a vector of \\p\\ free parameters
and the matrix is carried by a parameters7 parametrization whose free
values become ordinary scalar parameters, so every generic of the
`distrib` class applies to the object, with `theta` a flat named list of
numbers. Both default to an unstructured matrix in the log-Cholesky
parametrization, which has `p * (p + 1) / 2` free values.

## Usage

``` r
mvgaussian1_distrib(n_dim, sigma = NULL)

mvgaussian2_distrib(n_dim, omega = NULL)
```

## Arguments

- n_dim:

  The dimension \\p\\ of one observation. A single positive whole
  number, finite and at least 1. Any other value signals an error.

- sigma:

  A parameters7 parametrization of the covariance, of dimension `n_dim`
  and of full rank. Defaults to `parameters7::log_cholesky(n_dim)`.

- omega:

  A parameters7 parametrization of the precision, of dimension `n_dim`
  and of full rank. Defaults to `parameters7::log_cholesky(n_dim)`.

## Value

An S7 object of class
[MvGaussian1Distrib](https://statmodels7.github.io/distributions7/reference/MvGaussian1Distrib.md)
or
[MvGaussian2Distrib](https://statmodels7.github.io/distributions7/reference/MvGaussian1Distrib.md),
whose `param` is the parametrization supplied and whose `inverted`
records which side it carries. Its `params` are `mu1`, ..., `mup`
followed by the prefixed free names, `n_params` is `p + param@n_free`,
every entry of `params_bounds` is \\(-\infty, \infty)\\ and every link
is the identity.

## The two families

The numbering is the package's convention, the one that separates
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
and
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md):
one name per parametrization, because in a modeling framework a
parametrization determines what a linear predictor acts on. Here the
choice can change the model as well as the coordinates. A matrix
parametrization is closed under inversion when the inverse of every
matrix that it produces is one of its own;
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html),
[`parameters7::matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.html),
[`parameters7::diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.html),
[`parameters7::scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.html),
[`parameters7::compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.html)
and
[`parameters7::dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.html)
are, and
[`parameters7::correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.html),
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)
and
[`parameters7::autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.html)
are not.

Where it is closed the two families describe the same laws in different
coordinates, and a fit reaches the same maximum. Where it is not closed
they are different models: the inverse of an AR(1) covariance is
tridiagonal and is not an AR(1) matrix at any parameters, so
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)
on the covariance and on the precision give different maximized
log-likelihoods.

Both constructors accept any parametrization of the right dimension and
of full rank. `mvgaussian2_distrib(p, parameters7::ar1(p))` is the model
whose precision has the entries \\\sigma^2\rho^{\|i-j\|}\\, with
\\\sigma^2\\ the scale of `ar1()`, and it is a valid law. The AR(1)
process written on the precision side is
`mvgaussian2_distrib(p, parameters7::ar1_inv(p))`, whose parametrization
is that inverse.

## The arithmetic of the two sides

The derivatives of both families are computed from the derivative arrays
of \\\Sigma\\ and of \\\Sigma^{-1}\\. The parametrization supplies one
of the two, and
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md)
obtains the other: through
[`parameters7::param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.html)
and its siblings under a covariance parametrization, and by the chain
rule for an inverse under a precision one.

## Parameters

The mean contributes `mu1`, ..., `mup`. The matrix contributes the
parametrization's free values, prefixed by the matrix that they
describe: `sigma_` for a covariance and `omega_` for a precision. A
two-dimensional gaussian on an unstructured covariance therefore has the
five parameters `mu1`, `mu2`, `sigma_log_L1`, `sigma_log_L2` and
`sigma_L2.1`, while the same parametrization on the precision gives
`omega_log_L1` and the rest. The prefix distinguishes the two models in
a printed table, since the name of a free value identifies the
construction of the matrix and does not identify which matrix it is.

All of these parameters are unconstrained, so all of the links are the
identity. The constraint that makes the matrix positive definite is
carried by the matrix parametrization, which needs no link to express
it. The parameter scale and the link scale therefore coincide, and
`scale = "link"` returns the same numbers as `scale = "parameter"`.

## Reading a fit

The free values are coordinates of the parametrization and have no
direct interpretation.
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
carries the variance matrix of a fit onto the standard deviations and
the correlations by the delta method, and
[`print()`](https://rdrr.io/r/base/print.html) shows those. A precision
parametrization reports instead the conditional variances and the
partial correlations, which it describes directly: \\1/\Omega\_{jj} =
\operatorname{Var}(Y_j \mid Y\_{-j})\\, whose ratio to the marginal
variance is \\1 - R_j^2\\ for the regression of that coordinate on all
the others.

## Rank

A rank-deficient parametrization is rejected with an error. A singular
covariance gives a law supported on a subspace, with no density against
Lebesgue measure; a singular precision gives a quadratic form that is
flat along its null space and does not normalize. Such a parametrization
can serve as a penalty but does not define a density.

## The response

The response is an \\n \times p\\ matrix with one row per observation. A
plain vector of length \\p\\ is read as a single observation. A
parameter may not vary from observation to observation here, because the
parametrization describes one matrix for the whole sample.

## The distribution

\$\$f(y) = (2\pi)^{-p/2}\lvert \Sigma
\rvert^{-1/2}\exp\\\left\\-\tfrac{1}{2}(y-\mu)^\top\Sigma^{-1}(y-\mu)\right\\\$\$
on \\y \in \mathbb{R}^{p}\\, with \$\$\mathbb{E}\[Y\] = \mu, \qquad
\operatorname{Var}(Y) = \Sigma.\$\$

## Notation

\\\mu\\ is the mean vector, \\\Sigma\\ the covariance, \\\Omega\\ the
precision, \\p\\ the dimension of one observation and \\n\\ the number
of observations. \\\eta\\ is the free vector of the matrix
parametrization, the unconstrained coordinates that an optimizer moves.

## See also

[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the one-dimensional case,
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
for the heavy-tailed sibling,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the quantities that a fit reports,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate one,
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
for the default parametrization, and
[`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html)
for writing a structure on the other side.

## Examples

``` r
d <- mvgaussian1_distrib(2)
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

# The covariance that the free values describe is positive definite for any
# five numbers.
S <- mv_sigma(d, theta)
round(S, 4)
#>        v1     v2
#> v1 1.2214 0.4421
#> v2 0.4421 0.8303
eigen(S, only.values = TRUE)$values
#> [1] 1.5092462 0.5424766

# The density is compared with the formula written out.
y <- rbind(c(0, 0), c(1, -1))
mu <- c(0.5, -0.3)
ref <- apply(y, 1, function(r)
  -0.5 * (2 * log(2 * pi) + log(det(S)) + t(r - mu) %*% solve(S) %*% (r - mu)))
all.equal(distrib_pdf(d, y, theta, log = TRUE), as.numeric(ref))
#> [1] TRUE

# A structured covariance spends fewer parameters: two variances, no
# correlation.
mvgaussian1_distrib(2, parameters7::diagonal_matrix(2))@params
#> [1] "mu1"          "mu2"          "sigma_log_d1" "sigma_log_d2"

# The precision form is the other family, and its free values are prefixed
# by the matrix that they describe.
o <- mvgaussian2_distrib(2)
o@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
             omega_log_L2 = -0.2, omega_L2.1 = 0.4)
Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
all.equal(mv_sigma(o, th_o), solve(Om), check.attributes = FALSE)
#> [1] TRUE

# The conditional variance of the first coordinate given the second appears
# beside the marginal variance.
c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#> conditional    marginal 
#>   0.8187308   1.0141552 

# An AR(1) covariance and an AR(1) precision are different models, and the
# name of the family records which one is meant.
c(covariance = mvgaussian1_distrib(4, parameters7::ar1(4))@distrib_name,
  precision = mvgaussian2_distrib(4, parameters7::ar1(4))@distrib_name)
#>                              covariance                               precision 
#> "multivariate gaussian [4d, sigma=ar1]" "multivariate gaussian [4d, omega=ar1]" 
```
