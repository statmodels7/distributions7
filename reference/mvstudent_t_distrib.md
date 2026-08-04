# Construct a Multivariate Student's t Distribution

The elliptical \\t\\ on \\\mathbb{R}^p\\, with density \$\$f(y) \propto
\|\Sigma\|^{-1/2} \left(1 + \frac{(y-\mu)^\top \Sigma^{-1}
(y-\mu)}{\nu}\right)^{-(\nu+p)/2}.\$\$

## Usage

``` r
mvstudent_t_distrib(n_dim, sigma = NULL, link_nu = linkfunctions7::log_link())
```

## Arguments

- n_dim:

  The dimension \\p\\.

- sigma:

  A parameters7 structure for the scale matrix. Defaults to
  `parameters7::log_cholesky(n_dim)`.

- link_nu:

  The link for the degrees of freedom. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An object of class
[`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md).

## Details

\\\Sigma\\ is the **scale** matrix and not the covariance: the
covariance is \\\nu\Sigma/(\nu-2)\\ where it exists at all, and for
\\\nu \le 2\\ it does not while the distribution is perfectly well
defined.
[`mv_sigma`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
returns the scale matrix, which is the thing the parametrisation
carries, and
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md)
returns the covariance, which is a moment. Keeping the two apart is what
lets a fit run at \\\nu = 1.5\\.

**Parameters.** The mean contributes `mu1`, ..., `mup`, the structure
contributes its free values under the `sigma_` prefix, and `nu` is added
last. The mean and the matrix parameter are unconstrained and carry
identity links; `nu` is positive and carries a log link by default, so
unlike the multivariate gaussian this family's link scale is not its
parameter scale.

**Reading a fit.**
[`mv_summary`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
reports the square roots of the diagonal of the **scale** matrix and the
correlations. The correlations are the response's as well, a positive
multiple of a matrix leaving them alone, but the diagonal quantities are
not standard deviations of the response and are named `scale_sd_` to say
so.

**What it is for.** A gaussian fitted to data with a few outlying rows
inflates its covariance to cover them. A \\t\\ with \\\nu\\ estimated
does not: the observations far from the centre get a weight
\\(\nu+p)/(\nu+q)\\ that falls away with their Mahalanobis distance
\\q\\, which is what appears in the score below and what makes the fit
resistant. The gaussian is the limit \\\nu \to \infty\\.

**The expected information** has no closed form here and is approximated
by sampling.
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
therefore accepts `approx`, which it refuses for a family that computes
it exactly.

## See also

[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md),
[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)

## Examples

``` r
d <- mvstudent_t_distrib(2)
d@params
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
#> [6] "nu"          

theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4, nu = 5)

# the scale matrix is what the matrix parameter carries; the covariance is a moment
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.40
#> v2 0.4 1.16
variance(d, theta)
#>           v1        v2
#> v1 1.6666667 0.6666667
#> v2 0.6666667 1.9333333

# and below two degrees of freedom the covariance does not exist at all,
# while the density does
theta$nu <- 1.5
distrib_pdf(d, c(0, 0), theta)
#> [1] 0.1591549
variance(d, theta)
#>     v1  v2
#> v1 Inf Inf
#> v2 Inf Inf
```
