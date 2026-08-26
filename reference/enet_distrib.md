# Elastic-Net Distribution Object

Builds the density whose negative logarithm is the elastic-net penalty:
a Laplace and a Gaussian at the same location, multiplied together and
normalized. It exists so that `penalties7::elasticnet_penalty()` is the
same construction as ridge and lasso, a `penalties7::distrib_penalty()`
over a
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
family, instead of a branch of its own.

The parametrization is the one a reader of the penalty expects: an
overall rate \\\lambda\\ and a mixing weight \\\alpha\\, with \\a =
\lambda\alpha\\ the Laplace rate and \\c = \lambda(1-\alpha)\\ the
Gaussian one.

## Usage

``` r
enet_distrib(
  link_mu = identity_link(),
  link_lambda = log_link(),
  link_alpha = logit_link()
)
```

## Arguments

- link_mu:

  A `linkfunctions7` link object for the location \\\mu\\, which is
  unconstrained. Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_lambda:

  A link object for the overall rate \\\lambda\\, which must be strictly
  positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_alpha:

  A link object for the mixing weight \\\alpha\\, which must lie
  strictly inside \\(0, 1)\\. Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html).

## Value

An S7 object of class
[EnetDistrib](https://statmodels7.github.io/distributions7/reference/EnetDistrib.md),
inheriting from `continuous_distrib`. Its `params` are
`c("mu", "lambda", "alpha")`, its `bounds` `c(-Inf, Inf)`, its
`params_smooth` `c(mu = FALSE, lambda = TRUE, alpha = TRUE)`, and its
`link_params` the three links given here.

## The normalizing constant

\\Z = 2M(a/\sqrt{c})/\sqrt{c}\\ with \\M\\ the Mills ratio, finite at
both ends: \\2/a\\ as \\\alpha \to 1\\ and \\\sqrt{2\pi/c}\\ as \\\alpha
\to 0\\. Every derivative in the two rates is a polynomial in \\x =
a/\sqrt c\\ and \\G = \mathrm{d}\log M/\mathrm{d}x\\, with \\G' = 1 +
xG - G^2\\ closing the recursion, and the chain to \\(\lambda, \alpha)\\
is bilinear.

## What is closed form

All four derivative orders in the parameters, the expected information,
the distribution function, the quantile function, the mean, the variance
and the skewness. The **kurtosis** is not, and comes from the base class
by quadrature.

The expected information is closed form for a reason worth naming: in
the two rates the density is an exponential family with sufficient
statistics \\-\|z\|\\ and \\-z^2/2\\, so the information there is the
Hessian of \\\log Z\\, which the observed Hessian already computes. See
[`distrib_expected_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.EnetDistrib.md).

## The kink

Like the Laplace, the log-likelihood is not differentiable in \\\mu\\.
`params_smooth` records `mu = FALSE`, the observed curvature there is
\\-c\\ and misses the point mass at the location, and the information is
defined as the variance of the score instead. The two agree at every
other entry.

## Two numerical cautions

Both were found by sweeping the parameters, and both bite at ordinary
settings.

- \\\log M(x)\\ written directly adds \\x^2/2\\ to a log-probability of
  the same size and opposite sign, losing a digit per factor of ten in
  \\x\\. At \\\alpha = 1 - 10^{-12}\\ the argument reaches \\10^6\\ and
  the density was wrong in the fourth digit; past \\x = 30\\ the
  asymptotic series of the Mills ratio is used instead.

- \\\Phi(-x)\\ is exactly zero past \\x = 38\\, and the distribution
  function and the quantile divide by it. At \\\lambda = 20\\, \\\alpha
  = 0.995\\ the argument is 63, so both work through `log.p` and
  `qnorm(..., log.p = TRUE)`.

## Comparison with glmnet

`glmnet` standardizes the response internally for the Gaussian family
and corrects `lambda` by a single factor. That correction is exact for
the \\\ell_1\\ part, homogeneous of degree one, and cannot be for the
\\\ell_2\\ part, homogeneous of degree two. So a nominal \\(\lambda,
\alpha)\\ names two different objectives in the two packages whenever
\\\alpha \< 1\\. At \\\alpha = 1\\ the two agree to \\1.4\times10^{-4}\\
whatever the scale of the response.

## Parameter domains

- \\\mu \in (-\infty, \infty)\\

- \\\lambda \in (0, \infty)\\

- \\\alpha \in (0, 1)\\, open at both ends: at 1 it is
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
  and at 0
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
  both of which remain families of their own.

## The distribution

\$\$f(y) = \frac{1}{Z}\exp\\\left\\-a\|y-\mu\| -
\frac{c}{2}(y-\mu)^{2}\right\\, \qquad a = \lambda\alpha, \quad c =
\lambda(1-\alpha), \quad Z = \frac{2M(a/\sqrt{c})}{\sqrt{c}}\$\$ on \\y
\in \mathbb{R}\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \frac{1 +
xG}{c}, \qquad x = \frac{a}{\sqrt{c}}, \quad G = \frac{\mathrm{d}\log
M}{\mathrm{d}x}\$\$

## References

Zou, H. and Hastie, T. (2005). Regularization and variable selection via
the elastic net. *Journal of the Royal Statistical Society, Series B*
67, 301-320.

## See also

[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the two ends, `penalties7::elasticnet_penalty()` for the consumer,
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
for holding the location at zero, which is how a penalty uses this, and
[EnetDistrib](https://statmodels7.github.io/distributions7/reference/EnetDistrib.md)
for the class and its method list.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 0, lambda = 2, alpha = 0.5)

distrib_pdf(d, c(-1, 0, 1), th)
#> [1] 0.1701518 0.7625676 0.1701518

# It integrates to one, and the two ends are the two families.
integrate(function(u) distrib_pdf(d, u, th), -Inf, Inf)$value
#> [1] 1
rbind(enet = c(distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
                                        alpha = 1 - 1e-10)),
               distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
                                        alpha = 1e-10))),
      limit = c(distrib_pdf(laplace2_distrib(), 0.7,
                            list(mu = 0, lambda = 2)),
                dnorm(0.7, 0, 1 / sqrt(2))))
#>           [,1]      [,2]
#> enet  0.246597 0.3456374
#> limit 0.246597 0.3456374

# The kurtosis sits between the two ends, and is the one moment that is
# not closed form.
c(gaussian = 0, enet = kurtosis(d, th), laplace = 3)
#>  gaussian      enet   laplace 
#> 0.0000000 0.7658603 3.0000000 

# The location carries a kink, so its observed curvature and its
# information differ; every other entry agrees.
rbind(observed = unlist(distrib_hessian(d, 0, th)),
      expected = unlist(distrib_expected_hessian(d, 0, th)))
#>              mu_mu lambda_lambda alpha_alpha mu_lambda mu_alpha lambda_alpha
#> observed -1.000000    -0.1702646  -0.1159321         0        0   0.24452821
#> expected -2.525135    -0.1702646  -0.1159321         0        0  -0.04317471

# Held at zero, this is the prior a penalty is written from.
p <- fixed(d, mu = 0)
p@params
#> [1] "lambda" "alpha" 
```
