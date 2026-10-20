# Construct a Multivariate Student t Distribution

Two families on \\\mathbb{R}^p\\, differing in which matrix their
parametrization describes: `mvstudent_t1_distrib()` parametrizes the
scale matrix \\\Sigma\\ and `mvstudent_t2_distrib()` its inverse
\\\Sigma^{-1}\\. Each carries a location vector of \\p\\ free
parameters, the matrix through a parameters7 parametrization, and a
positive degrees-of-freedom parameter. The free values become ordinary
scalar parameters, so the generics of the `distrib` class accept the
object, with `theta` a flat named list of numbers; the distribution
function and the quantile function signal an error. Both default to an
unstructured matrix in the log-Cholesky parametrization, which has
`p * (p + 1) / 2` free values, and to a log link on `nu`.

## Usage

``` r
mvstudent_t1_distrib(n_dim, sigma = NULL, link_nu = linkfunctions7::log_link())

mvstudent_t2_distrib(n_dim, omega = NULL, link_nu = linkfunctions7::log_link())
```

## Arguments

- n_dim:

  The dimension \\p\\ of one observation. The value must be a single
  positive whole number, finite and at least 1. Any other value signals
  an error.

- sigma:

  A parameters7 parametrization of the scale matrix, of dimension
  `n_dim` and of full rank. Defaults to
  `parameters7::log_cholesky(n_dim)`. A rank-deficient parametrization
  is rejected because the density does not normalize.

- link_nu:

  The link carrying the degrees of freedom to the unconstrained scale, a
  [`linkfunctions7::link`](https://statmodels7.github.io/linkfunctions7/reference/link.html)
  object. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which keeps \\\nu\\ strictly positive at every point of the free
  scale.

- omega:

  A parameters7 parametrization of the inverse scale matrix, of
  dimension `n_dim` and of full rank. Defaults to
  `parameters7::log_cholesky(n_dim)`. It is the precision of the
  response only up to the factor \\(\nu-2)/\nu\\; see Details.

## Value

An S7 object of class
[MvStudentT1Distrib](https://statmodels7.github.io/distributions7/reference/MvStudentT1Distrib.md)
or
[MvStudentT2Distrib](https://statmodels7.github.io/distributions7/reference/MvStudentT1Distrib.md),
whose `param` is the parametrization supplied and whose `inverted`
records which matrix `param` carries. Its `params` are `mu1`, ...,
`mup`, the prefixed free names, then `nu`; `n_params` is
`p + param@n_free + 1`; every `params_bounds` entry is \\(-\infty,
\infty)\\ except that of `nu`, which is \\(0, \infty)\\; and every link
is the identity except that of `nu`.

## The two matrices

Each of the two matrices is neither a covariance nor a precision.
\\\Sigma\\ is the scale matrix, the matrix in the quadratic form of the
density, and the moments are obtained from it by scaling in \\\nu\\:

\$\$\operatorname{Var}(Y) = \frac{\nu}{\nu-2}\\\Sigma \quad (\nu \> 2),
\qquad \operatorname{Var}(Y)^{-1} = \frac{\nu-2}{\nu}\\\Sigma^{-1}.\$\$

So the matrix of `mvstudent_t2_distrib()` equals the precision of the
response only up to the factor \\(\nu-2)/\nu\\, which is 2/3 at \\\nu =
6\\; a coordinate named `omega_log_L1` is a free value of
\\\Sigma^{-1}\\ and not of the precision. The prefixes follow the
convention of
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and name which of the two matrices a free value builds.

For \\\nu \le 2\\ the covariance does not exist while the density is
well defined, and for \\\nu \le 1\\ the mean does not exist either. For
that reason
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
returns the scale matrix, which the parametrization carries and which
exists at every admissible \\\nu\\, and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
returns the covariance, a moment, which is `Inf` in every entry for
\\\nu \le 2\\; [`mean()`](https://rdrr.io/r/base/mean.html) returns
`NaN` for \\\nu \le 1\\. A fit therefore runs at \\\nu = 1.5\\.

Two readings carry over from the matrix to the response without a
factor, because a positive multiple cancels out of a ratio:

- the correlations of \\\Sigma\\ are the correlations of \\Y\\;

- the partial correlations read off \\\Sigma^{-1}\\, as
  \\-\Omega\_{jk}/\sqrt{\Omega\_{jj}\Omega\_{kk}}\\, are the partial
  correlations of \\Y\\.

The conditional variance does not carry over. For a gaussian
\\1/\Omega\_{jj}\\ is the conditional variance of \\Y_j\\ given the
rest; here it is the Schur complement of the scale matrix, so it is a
conditional scale, and the conditional law of a \\t\\ is itself a \\t\\
whose scale carries a further factor depending on the conditioning
values.
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
names that quantity `cscale_vj` and not `cvar_vj` for this reason, as it
names the diagonal quantities `scale_sd_vj`.

## Which side to write a structure on

Where the matrix parametrization is closed under inversion the two
families describe the same laws in different coordinates. Where it is
not they are different models:
`mvstudent_t2_distrib(p, parameters7::ar1(p))` is the family whose
inverse scale matrix has the entries \\\sigma^2\rho^{\|i-j\|}\\, with
\\\sigma^2\\ the scale of `ar1()`, and it is a valid law. The AR(1)
pattern on the scale matrix, in the family that parametrizes the
inverse, is `mvstudent_t2_distrib(p, parameters7::ar1_inv(p))`, whose
parametrization is that inverse.

## Parameters

The location contributes `mu1`, ..., `mup`, the parametrization
contributes its free values under the `sigma_` prefix for the scale
matrix and the `omega_` prefix for its inverse, and `nu` comes last. The
location and the matrix parameters are unconstrained and carry identity
links; `nu` is positive and carries `link_nu`. Unlike
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md),
this family therefore has a link scale that differs from its parameter
scale: under the default log link, `scale = "link"` applies the chain
rule of the log link to every derivative in `nu`.

## The weight of an observation

When a gaussian is fitted to data with a few outlying rows, the
estimated covariance grows until it covers them. A \\t\\ whose \\\nu\\
is estimated does not: an observation at squared Mahalanobis distance
\\q\\ enters every derivative through the weight \\c =
(\nu+p)/(\nu+q)\\, which falls as \\1/q\\. The gaussian is the limit
\\\nu \to \infty\\, where \\c \equiv 1\\.

## The expected information

The expected information has a closed form, obtained from the family's
own scale mixture; see
[`distrib_expected_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md)
for the four blocks.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
therefore rejects an `approx` argument on
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for this family, as it does for any family that computes its expected
information exactly.

## Reading a fit

For the scale parametrization,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
reports the square roots of the diagonal of the scale matrix, named
`scale_sd_v1`, ..., `scale_sd_vp` because they are not standard
deviations of the response, and the correlations, which are the
response's as well. For the inverse parametrization it reports instead
the conditional scales `cscale_v1`, ..., `cscale_vp` and the partial
correlations `pcor_vj_vk`, which are the response's exactly.

## Degrees of freedom far from the gaussian limit and close to it

Every quantity that depends only on \\\nu\\ (the constant of the density
and its derivatives, and the expected information in \\\nu\\) is
computed by
[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
and
[`mvt_info_nunu()`](https://statmodels7.github.io/distributions7/reference/mvt_info_nunu.md)
without the cancellation of the direct difference of gamma and polygamma
functions, and the part that carries the data by
[`mvt_h_d()`](https://statmodels7.github.io/distributions7/reference/mvt_h_d.md);
all of them vanish as \\\nu\\ grows, while the terms of the direct forms
do not.

## The response

The response is an \\n \times p\\ matrix with one row per observation. A
plain vector of length \\p\\ is read as a single observation. A
parameter may not vary from observation to observation.

## The distribution

\$\$f(y) =
\frac{\Gamma\\\left(\frac{\nu+p}{2}\right)}{\Gamma\\\left(\frac{\nu}{2}\right)(\nu\pi)^{p/2}\lvert
\Sigma \rvert^{1/2}}\left(1 + \frac{q}{\nu}\right)^{-(\nu+p)/2}, \qquad
q = (y-\mu)^\top\Sigma^{-1}(y-\mu)\$\$ on \\y \in \mathbb{R}^{p}\\, with
\$\$\mathbb{E}\[Y\] = \mu \\ (\nu \> 1), \qquad \operatorname{Var}(Y) =
\frac{\nu}{\nu-2}\Sigma \\ (\nu \> 2).\$\$

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\p\\ the dimension of one observation and \\q\\ the
squared Mahalanobis distance of an observation from the location.
\\\eta\\ is the free vector of the matrix parametrization.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
for the limiting family,
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the one-dimensional case,
[`mv_marginal.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvStudentTDistrib.md)
for a marginal, which keeps the same \\\nu\\,
[`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html)
for writing a structure on the other side, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate one.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
d
#> Distribution: Multivariate Student T [2d, Sigma=log_cholesky]
#> Type:         Continuous, 2-dimensional
#> Dimensions:   multivariate
#> 
#> Parameters:
#>   mu1          (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   mu2          (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_log_L1 (scale)              | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_log_L2 (scale)              | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_L2.1   (scale)              | Link: identity   | Domain: (-Inf, Inf)
#>   nu           (degrees of freedom) | Link: log        | Domain: (0, Inf)

theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

# The parametrization carries the scale matrix; the covariance is a moment
# and is larger by nu / (nu - 2).
mv_sigma(d, theta)
#>           v1        v2
#> v1 1.2214028 0.4420684
#> v2 0.4420684 0.8303200
all.equal(variance(d, theta), (6 / 4) * mv_sigma(d, theta))
#> [1] TRUE

# For nu <= 2 the covariance is infinite and the density is finite; for
# nu <= 1 the mean does not exist.
theta$nu <- 1.5
c(density = distrib_pdf(d, c(0, 0), theta), variance = variance(d, theta)[1, 1])
#>   density  variance 
#> 0.1018402       Inf 
theta$nu <- 0.8
c(density = distrib_pdf(d, c(0, 0), theta), mean = mean(d, theta)[1])
#>    density    mean.v1 
#> 0.08456772        NaN 

# At four units of scale from the location, a t at nu = 6 has a far larger
# density than a gaussian with the same matrix.
theta$nu <- 6
g <- mvgaussian1_distrib(2)
c(t = distrib_pdf(d, c(4.5, -4.3), theta),
  gaussian = distrib_pdf(g, c(4.5, -4.3), theta[1:5]))
#>            t     gaussian 
#> 1.413233e-05 6.118878e-14 

# The other parametrization carries the inverse scale matrix, and the prefix
# of its free values names it.
o <- mvstudent_t2_distrib(2)
o@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
#> [6] "nu"          

# Its matrix is the precision of the response only up to (nu - 2) / nu.
th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
             omega_log_L2 = -0.2, omega_L2.1 = 0.4, nu = 6)
Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
all.equal(unname(solve(variance(o, th_o))),
          unname((6 - 2) / 6 * unclass(Om)))
#> [1] TRUE

# The correlations are the response's exactly: a positive multiple cancels
# out of a ratio.
all.equal(cov2cor(mv_sigma(o, th_o)), cov2cor(variance(o, th_o)))
#> [1] TRUE
```
