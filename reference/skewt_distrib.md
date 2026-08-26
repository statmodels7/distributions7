# Skew t Distribution Object

Builds a skew \\t\\ distribution object with location \\\mu\\, scale
\\\sigma\\, shape \\\alpha\\ and degrees of freedom \\\nu\\. It is the
four-parameter family a location-scale-shape framework wants: the scale,
the asymmetry and the tail weight are three separate parameters, each of
which can carry its own linear predictor.

Three families sit inside it. \\\alpha = 0\\ is
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
large \\\nu\\ approaches
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
and both together give the Gaussian.

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

  A `linkfunctions7` link object for the location \\\mu\\, which is
  unconstrained. Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link object for the scale \\\sigma\\, which must be strictly
  positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_alpha:

  A link object for the shape \\\alpha\\, which is unconstrained.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_nu:

  A link object for the degrees of freedom \\\nu\\, which must be
  strictly positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class
[SkewTDistrib](https://statmodels7.github.io/distributions7/reference/SkewTDistrib.md),
inheriting from `continuous_distrib`. Its `params` are
`c("mu", "sigma", "alpha", "nu")`, its `bounds` `c(-Inf, Inf)`, and its
`link_params` the four links given here.

## Why the fourth parameter

The skew normal of
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
reaches a skewness of at most 0.99527 and an excess kurtosis of at most
0.86918, both approached only as \\\|\alpha\| \to \infty\\. Adding
\\\nu\\ removes both ceilings: measured at \\\alpha = 50\\, the skewness
is 1.190 at \\\nu = 30\\, 2.050 at \\\nu = 6\\ and 3.998 at \\\nu = 4\\.

## What is closed form and what is not

Every derivative in \\(\mu, \sigma, \alpha)\\ is closed form. Everything
involving \\\nu\\ is not, because the density contains \\T\_{\nu+1}\\
and the derivative of a Student \\t\\ distribution function with respect
to its degrees of freedom has no elementary expression. It is the same
obstruction that stops the gamma and beta distribution functions from
having closed-form shape derivatives.

Those components come from a single stencil applied to an analytic
quantity, never from a difference of a difference:

|  |  |  |
|----|----|----|
| **component** | **route** | **agreement with an independent route** |
| \\\mu, \sigma, \alpha\\ (score) | closed form | \\10^{-12}\\ |
| \\\nu\\ (score) | [`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md) on \\\ell\\ | \\5\times10^{-11}\\ |
| \\(\mu,\sigma,\alpha)\\ block (Hessian) | closed form | \\10^{-12}\\ |
| \\\nu\\ with another parameter | [`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md) on the analytic score | to the printed digit |
| \\\nu\\ twice | [`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md) on \\\ell\\ | \\2\times10^{-9}\\ |
| \\\nu\\ three times | [`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md) on \\\ell\\ | \\10^{-4}\\ |
| \\\nu\\ four times | [`fd5_fourth()`](https://statmodels7.github.io/distributions7/reference/fd5_fourth.md) on \\\ell\\ | about one figure |

## The tolerance a fit can ask for

The score in \\\nu\\ cannot be computed more accurately than that table,
so no stopping rule on the gradient can be satisfied below it however
good the optimizer is.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
takes its rule from the method it is given, and `crit_grad()`'s default
tolerance of \\10^{-6}\\ is tested on the score **per observation**,
which leaves room: measured on samples of 500 to 4000, the run converges
with a summed score between \\1.5\times10^{-7}\\ and
\\2.4\times10^{-5}\\, i.e. between \\3\times10^{-10}\\ and
\\9\times10^{-9}\\ per observation.

## Fitting

The expected information has no closed form, so it is approximated by
the strategy named in `approx`, at one quadrature per component.
`method = optimizers7::newton()` is much the cheaper route: the observed
Hessian is the closed form above and needs no integration.

The distribution function and the quantile function likewise have no
elementary form; the base class integrates the density and inverts the
result by root finding.

## Moments

They exist only up to order \\\nu\\: the mean requires \\\nu \> 1\\, the
variance \\\nu \> 2\\, the skewness \\\nu \> 3\\ and the excess kurtosis
\\\nu \> 4\\, and each returns `NaN` below its threshold. The density is
well defined at every positive \\\nu\\, which is why the moments and the
parameters are kept apart.

## Identification near symmetry

The skew normal's expected information loses a rank at \\\alpha = 0\\,
because its shape score is a fixed multiple of its location score there.
**The skew \\t\\ does not inherit that at finite \\\nu.\\** The tilting
argument carries \\c = \sqrt{(\nu+1)/(\nu+z^2)}\\, which depends on the
observation, so the two scores are not proportional: measured at
\\\alpha = 0\\, \\\nu = 6\\, their ratio spreads over 0.76 across
observations, where the skew normal's is constant to the last digit. The
spread falls as \\O(1/\nu)\\, so the singularity is inherited only in
the limit: 0.068 at \\\nu = 100\\, \\7\times10^{-4}\\ at \\10^4\\.

What is weakly identified here is \\\nu\\ itself. The smallest
eigenvalue of the information belongs to that direction and falls with
\\\nu\\: measured at \\\alpha = 0\\, it is \\3.6\times10^{-3}\\ at \\\nu
= 4\\, \\8.4\times10^{-4}\\ at 6 and \\1.7\times10^{-6}\\ at 30. A
nearly Gaussian tail carries little information about how heavy it is.

## Parameter domains

- \\\mu \in (-\infty, \infty)\\

- \\\sigma \in (0, \infty)\\

- \\\alpha \in (-\infty, \infty)\\

- \\\nu \in (0, \infty)\\

## References

Azzalini, A. and Capitanio, A. (2003). Distributions generated by
perturbation of symmetry with emphasis on a multivariate skew t
distribution. *Journal of the Royal Statistical Society, Series B* 65,
367-389.

## See also

[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the \\\nu \to \infty\\ limit,
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the \\\alpha = 0\\ case,
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
for the scalar functions the derivatives are built from, and
[SkewTDistrib](https://statmodels7.github.io/distributions7/reference/SkewTDistrib.md)
for the class and its method list.

## Examples

``` r
d <- skewt_distrib()
d@params
#> [1] "mu"    "sigma" "alpha" "nu"   
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 5)

distrib_pdf(d, c(-1, 0, 1), th)
#> [1] 0.005274116 0.379606690 0.434085479

# Shape zero is the Student t.
all.equal(distrib_pdf(d, c(-1, 0, 1),
                      list(mu = 0, sigma = 1, alpha = 0, nu = 5)),
          dt(c(-1, 0, 1), df = 5))
#> [1] TRUE

# It passes both of the skew normal's ceilings.
rbind(skew_t = c(skewness(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5)),
                 kurtosis(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5))),
      skew_normal_bound = c(0.99527, 0.86918))
#>                       [,1]     [,2]
#> skew_t            2.477531 19.42035
#> skew_normal_bound 0.995270  0.86918

# Moments exist only up to order nu.
t(vapply(c(1.5, 2.5, 3.5, 4.5), function(v) {
  p <- list(mu = 0, sigma = 1, alpha = 3, nu = v)
  c(nu = v, mean = mean(d, p), var = variance(d, p),
    skew = skewness(d, p), kurt = kurtosis(d, p))
}, numeric(5)))
#>       nu      mean      var     skew     kurt
#> [1,] 1.5 1.9394975      NaN      NaN      NaN
#> [2,] 2.5 1.1441396 3.690944      NaN      NaN
#> [3,] 3.5 0.9875438 1.358091 6.021607      NaN
#> [4,] 4.5 0.9210146 0.951732 2.575662 34.78593

# The observed Hessian is the cheap route: this family has no closed-form
# expected information, so Fisher scoring would quadrature it at every step.
set.seed(1)
x <- distrib_rng(d, 200, th)
coef(fit_distrib(d, x, method = optimizers7::newton(), start = th))
#>          mu       sigma       alpha          nu 
#> -0.02273172  0.95605508  2.98317264  3.77234496 
```
