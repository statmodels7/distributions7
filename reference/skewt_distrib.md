# Skew t Distribution Object

Creates a distribution object for Azzalini's skew \\t\\ distribution,
with location \\\mu\\, scale \\\sigma\\, shape \\\alpha\\ and degrees of
freedom \\\nu\\. It contains the Student \\t\\ (\\\alpha = 0\\), the
skew normal (\\\nu \to \infty\\) and the gaussian (both).

## Usage

``` r
skewt_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_alpha = identity_link(),
  link_nu = log_link()
)
```

## Arguments

- link_mu:

  A link function object for the location \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link function object for the scale \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_alpha:

  A link function object for the shape \\\alpha\\, which is
  unconstrained. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_nu:

  A link function object for the degrees of freedom \\\nu\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class
[`SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/SkewTDistrib.md)
(inheriting from `continuous_distrib`).

## Details

This is the four-parameter family a location-scale-shape framework
wants: the scale, the skewness and the tail weight are three separate
parameters, each of which can be given its own linear predictor. The
skew normal of
[`skewnormal_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md)
can reach a skewness of at most \\0.995\\ and an excess kurtosis of at
most \\0.87\\; adding \\\nu\\ removes both bounds.

**Probability density function**, with \\z = (y-\mu)/\sigma\\ and \\w =
\alpha z\sqrt{(\nu+1)/(\nu+z^2)}\\: \$\$f(y; \mu, \sigma, \alpha, \nu) =
\dfrac{2}{\sigma}\\t\_\nu(z)\\T\_{\nu+1}(w)\$\$

**What is closed form and what is not.** The score and the observed
Hessian are closed form in \\(\mu, \sigma, \alpha)\\. Everything
involving \\\nu\\ is not, because the density contains \\T\_{\nu+1}\\,
whose derivative with respect to its degrees of freedom has no
elementary expression — the same obstruction that stops the gamma and
beta distribution functions from having closed-form shape derivatives.
Those components come from a single finite-difference stencil applied to
an analytic quantity, never from a difference of a difference:

|  |  |  |
|----|----|----|
| **component** | **route** | **error, summed over n** |
| \\\mu, \sigma, \alpha\\ (score) | closed form | machine precision |
| \\\nu\\ (score) | five-point stencil on \\\ell\\ | \\10^{-11}\\ to \\10^{-9}\\ |
| \\(\mu,\sigma,\alpha)\\ block (Hessian) | closed form | machine precision |
| \\\nu\\ with another parameter | five-point stencil on the analytic score | about \\10^{-8}\\ |
| \\\nu\\ twice | five-point stencil on \\\ell\\ | about \\10^{-6}\\ |

**The tolerance a fit can ask for.** The score in \\\nu\\ cannot be
computed more accurately than the table above, so a stopping rule on the
gradient cannot be satisfied below that level however good the optimiser
is.
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
tests the score **per observation**, and its default of \\10^{-10}\\
leaves room: on samples of 500 to 4000 the summed score reaches
\\10^{-10}\\ to \\3 \times 10^{-9}\\, which is \\10^{-13}\\ per
observation, and the fit converges in four or five iterations. A rule
expressed on the summed gradient at that tolerance would not be
attainable, which is why the tolerance is not expressed that way.

**Distribution function.** There is no elementary form, so the base
class integrates the density and inverts the result by root finding.

**Expected information** has no closed form either and is approximated
by the strategy named in `approx`, the default being the score variance.
That approximation costs one quadrature per component, so
`method = "newton"` is much the cheaper way to fit this family: the
observed Hessian is the closed form above and needs no integration.

**Moments** exist only up to order \\\nu\\: the mean requires \\\nu \>
1\\, the variance \\\nu \> 2\\, the skewness \\\nu \> 3\\ and the
kurtosis \\\nu \> 4\\, and each returns `NaN` below its threshold. The
density is perfectly well defined there, which is why the moments and
the parameters are kept apart.

**Special cases.** \\\alpha = 0\\ is
[`student_t_distrib`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md);
large \\\nu\\ approaches
[`skewnormal_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md).
The information is singular in \\\alpha\\ at \\\alpha = 0\\ for the same
reason as in the skew normal.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

- \\\alpha \in (-\infty, +\infty)\\

- \\\nu \in (0, +\infty)\\

## References

Azzalini, A. and Capitanio, A. (2003). Distributions generated by
perturbation of symmetry with emphasis on a multivariate skew t
distribution. *Journal of the Royal Statistical Society, Series B* 65,
367-389.

## Examples

``` r
d <- skewt_distrib()
d@params
#> [1] "mu"    "sigma" "alpha" "nu"   

theta <- list(mu = 0, sigma = 1, alpha = 3, nu = 5)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.005274116 0.379606690 0.434085479
distrib_gradient(d, c(-1, 0, 1), theta)
#> $mu
#> [1] -4.2263842 -2.5155765  0.9607996
#> 
#> $sigma
#> [1]  3.2263842 -1.0000000 -0.0392004
#> 
#> $alpha
#> [1] -1.29055368  0.00000000  0.01568016
#> 
#> $nu
#> [1] -0.185070253  0.009813847  0.021128296
#> 

# shape zero is the Student t
max(abs(distrib_pdf(d, c(-1, 0, 1), list(mu = 0, sigma = 1, alpha = 0, nu = 5)) -
        stats::dt(c(-1, 0, 1), df = 5)))
#> [1] 2.775558e-17

# the family reaches skewness the skew normal cannot
c(skew_t = skewness(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5)),
  skew_normal_bound = 0.9953)
#>            skew_t skew_normal_bound 
#>          2.477531          0.995300 

# The observed Hessian is the cheap route here: this family has no
# closed-form expected information, so Fisher scoring would approximate it
# by quadrature at every step.
set.seed(1)
y <- distrib_rng(d, 200, theta)
coef(fit_distrib(d, y, method = optimizers7::newton(), start = theta))
#>          mu       sigma       alpha          nu 
#> -0.02273189  0.95605544  2.98317484  3.77234947 
```
