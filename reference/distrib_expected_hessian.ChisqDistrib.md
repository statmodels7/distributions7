# Chi-Squared Expected Hessian

Returns the same number as
[`distrib_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md),
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\psi'(\mu/2)}{4},\$\$ the observed second derivative being free
of the response and so equal to its own expectation. Nothing is
averaged, integrated or simulated, so `approx` and `nsim` are ignored
and `y` is read only for its length. The value is negative at every
\\\mu\\, the trigamma function being positive, so the information is
positive throughout.

The identity holds on the **parameter** scale. On the link scale the
second-order chain rule adds \\h''(\eta)\\\partial\ell/\partial\mu\\ to
the observed Hessian; the expected version drops that term because the
score has mean zero, and a finite sample does not. Fisher scoring and
Newton's method therefore take different steps here and agree at the
optimum, where the summed score vanishes.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, the expectation being exact. Accepted so that the signature
  matches the generic's, where it selects between the Bartlett,
  quadrature, Monte Carlo and outer-product routes.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list with one numeric vector, `mu_mu`, of length
`max(length(y), length(mu))` and constant within itself when the
parameter is.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\mu^2\]\\, the expectation of the
**observed information** under the model. The chi-squared is a regular
family, so the second Bartlett identity holds and this equals the
variance of the score. \\\psi'\\ is the trigamma function.

## See also

[`distrib_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md),
which returns the same number;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step; and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
th <- list(mu = 4)

# A constant, and identical to the observed Hessian.
unique(distrib_expected_hessian(d, c(1, 4, 9), th)$mu_mu)
#> [1] -0.1612335
identical(distrib_expected_hessian(d, c(1, 4, 9), th),
          distrib_hessian(d, c(1, 4, 9), th))
#> [1] TRUE

# Negative at every mu, the trigamma function being positive.
vapply(c(0.5, 4, 40),
       function(m) distrib_expected_hessian(d, 0, list(mu = m))$mu_mu,
       numeric(1))
#> [1] -4.29933229 -0.16123352 -0.01281771

# The two fitting methods differ on the link scale and land together.
set.seed(7)
z <- distrib_rng(d, 2000, th)
c(newton = coef(fit_distrib(d, z, method = "newton")),
  fisher = coef(fit_distrib(d, z, method = "fisher")))
#> newton.mu fisher.mu 
#>  4.038376  4.038376 
```
