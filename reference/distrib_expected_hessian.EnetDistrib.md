# Elastic-Net Expected Information

Returns the expected Hessian in closed form. The family already carries
every piece of it: in the two rates the density is an exponential
family, so the information is the Hessian of its own log normalizing
constant, and `.enet_logz_derivs()` computes that for the observed
Hessian already.

`expected`, `approx` and `nsim` are all ignored, the answer being exact.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector. Its values do not enter the result, which is an
  expectation; only its length does, through recycling.

- theta:

  A named list with components `mu`, `lambda` and `alpha`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  Ignored: the answer is exact. Accepted so that the signature matches
  the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `mu_mu`, `lambda_lambda`, `alpha_alpha`, `mu_lambda`, `mu_alpha`,
`lambda_alpha`. The two mixed location entries are exactly zero.

## The rate block

With sufficient statistics \\-\|z\|\\ and \\-z^2/2\\, \\\log Z\\ is the
cumulant generating function and \$\$I\_{aa} =
\operatorname{Var}(\|z\|), \quad I\_{ac} =
\operatorname{Cov}\\\left(-\|z\|, -\tfrac{z^2}{2}\right), \quad I\_{cc}
= \tfrac{1}{4}\operatorname{Var}(z^2).\$\$ The map to \\(\lambda,
\alpha)\\ is bilinear, so the information transforms by \\J^\top I J\\
with no second-derivative term.

Two of the three rate entries of the observed Hessian carry no data and
are therefore their own expectations, which is why `lambda_lambda` and
`alpha_alpha` repeat them exactly. The third does not: `lambda_alpha`
carries \\-\|z\| + z^2/2\\, whose expectation \\\partial_a\log Z -
\partial_c\log Z\\ cancels the constant of the same value sitting beside
it and leaves the \\\log Z\\ term alone.

## The location, where the kink lives

The observed second derivative in \\\mu\\ is \\-c\\, which misses the
point mass \\\mathrm{d}\\\mathrm{sgn}(z)/\mathrm{d}z = 2\delta(z)\\ the
density carries at its own location, exactly as the Laplace does. The
information is defined there as the variance of the score, and
\$\$I\_{\mu\mu} = E\left\[(a\\\mathrm{sgn}(z) + cz)^2\right\] = a^2 +
2ac\\E\|z\| + c^2E\[z^2\] = a^2 - 2ac\\\partial_a \log Z - 2c^2
\partial_c \log Z,\$\$ since \\E\|z\| = -\partial_a\log Z\\ and
\\E\[z^2\] = -2\partial_c\log Z\\.

Measured at \\\lambda = 2\\: it is 3.99999996 as \\\alpha \to 1\\, which
is the Laplace's \\\lambda^2 = 1/\sigma^2\\, and exactly 2 as \\\alpha
\to 0\\, which is the Gaussian's \\c\\.

The location stays **orthogonal** to both rates, every cross-expectation
vanishing by symmetry, so `mu_lambda` and `mu_alpha` are exactly zero.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, \\z = y-\mu\\ and
\\Z\\ the normalizing constant.

## See also

[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md)
for the observed curvature, which differs from this only in the
location;
[`distrib_expected_hessian.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Laplace2Distrib.md)
for the same phenomenon in the family this one contains; and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 0, lambda = 2, alpha = 0.5)
e <- distrib_expected_hessian(d, 0, th)
names(e)
#> [1] "mu_mu"         "lambda_lambda" "alpha_alpha"   "mu_lambda"    
#> [5] "mu_alpha"      "lambda_alpha" 

# The location is orthogonal to both rates.
c(mu_lambda = e$mu_lambda, mu_alpha = e$mu_alpha)
#> mu_lambda  mu_alpha 
#>         0         0 

# The observed and expected curvature differ in the location and nowhere
# else among the entries that carry no data.
h <- distrib_hessian(d, 0, th)
rbind(observed = c(h$mu_mu, h$lambda_lambda, h$alpha_alpha),
      expected = c(e$mu_mu, e$lambda_lambda, e$alpha_alpha))
#>               [,1]       [,2]       [,3]
#> observed -1.000000 -0.1702646 -0.1159321
#> expected -2.525135 -0.1702646 -0.1159321

# The location entry reaches both ends: lambda^2 at alpha -> 1, c at 0.
inf_mu <- function(a)
  -distrib_expected_hessian(d, 0, list(mu = 0, lambda = 2, alpha = a))$mu_mu
rbind(alpha = c(1 - 1e-8, 0.5, 1e-8),
      information = vapply(c(1 - 1e-8, 0.5, 1e-8), inf_mu, 0),
      limit = c(2^2, NA, 2 * (1 - 1e-8)))
#>             [,1]     [,2]  [,3]
#> alpha          1 0.500000 1e-08
#> information    4 2.525135 2e+00
#> limit          4       NA 2e+00

# The strategy argument is ignored, the answer being exact.
identical(e, distrib_expected_hessian(d, 0, th, approx = "mc", nsim = 5))
#> [1] TRUE
```
