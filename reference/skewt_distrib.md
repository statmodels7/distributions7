# Skew t Distribution Object

Builds a skew \\t\\ distribution object with location \\\mu\\, scale
\\\sigma\\, shape \\\alpha\\ and degrees of freedom \\\nu\\. The scale,
the asymmetry and the tail weight are three separate parameters, and
each of them can carry its own linear predictor.

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

## Range of the skewness and the kurtosis

The skew normal of
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
reaches a skewness of at most 0.99527 and an excess kurtosis of at most
0.86918, both approached only as \\\|\alpha\| \to \infty\\. Adding
\\\nu\\ removes both ceilings: at \\\alpha = 50\\ the skewness is 1.12
at \\\nu = 30\\, 2.05 at \\\nu = 6\\ and 4.00 at \\\nu = 4\\, and it
grows without bound as \\\nu\\ decreases to 3.

## How the derivatives are computed

The density contains \\T\_{\nu+1}\\, and the derivatives of a Student
\\t\\ distribution function in its degrees of freedom have no elementary
expression. They are integrals of the derivatives of the \\t\\ density
and are taken by quadrature. The derivatives to order five are closed
form in \\\mu\\, \\\sigma\\ and \\\alpha\\, and those in \\\nu\\ come
from that quadrature; see
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md).

## Fitting

The expected information has no closed form and is computed by one
quadrature over \\z\\ per distinct \\(\alpha, \nu)\\, so Fisher scoring
on an intercept-only fit pays for one integral per iteration.

The distribution function and the quantile function likewise have no
elementary form. The distribution function integrates the density by the
compiled rule of
[`compiled_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md),
and the quantile function inverts it by root finding in
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## Moments

They exist only up to order \\\nu\\: the mean requires \\\nu \> 1\\, the
variance \\\nu \> 2\\, the skewness \\\nu \> 3\\ and the excess kurtosis
\\\nu \> 4\\, and each returns `NaN` below its threshold. The density is
well defined at every positive \\\nu\\, which is why the moments and the
parameters are kept apart.

## Identification near symmetry

The skew normal's expected information loses a rank at \\\alpha = 0\\,
because its shape score is a fixed multiple of its location score there.
For the skew \\t\\ at a finite \\\nu\\ the ratio of the score in
\\\alpha\\ to the score in \\\mu\\ at \\\alpha = 0\\ is proportional to
\\\sqrt{(\nu+z^2)/(\nu+1)}\\, which depends on the observation through
\\z\\, whereas the skew normal's ratio is a constant. The relative
variation of the ratio over \\z\\ is of order \\z^2/\nu\\, so the
information becomes singular only as \\\nu \to \infty\\.

The parameter that is weakly identified is \\\nu\\. At \\\alpha = 0\\
the smallest eigenvalue of the expected information corresponds to the
direction of \\\nu\\ and decreases as \\\nu\\ grows, because a tail
close to the Gaussian one carries little information about its weight.

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
#> [1,] 1.5 1.9394975      Inf      NaN      NaN
#> [2,] 2.5 1.1441396 3.690944      NaN      Inf
#> [3,] 3.5 0.9875438 1.358091 6.021607      Inf
#> [4,] 4.5 0.9210146 0.951732 2.575662 34.78593

# Fitting by Fisher scoring, the default, and by Newton on the observed
# Hessian reaches the same estimate.
set.seed(1)
x <- distrib_rng(d, 200, th)
rbind(fisher = coef(fit_distrib(d, x, start = th)),
      newton = coef(fit_distrib(d, x, method = optimizers7::newton(), start = th)))
#>                 mu     sigma    alpha       nu
#> fisher -0.02273143 0.9560544 2.983171 3.772340
#> newton -0.02273172 0.9560551 2.983173 3.772345
```
