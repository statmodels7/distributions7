# Zero-Adjusted Continuous Expected Information

Computes the expectation of the observed Hessian in closed form:
\$\$\mathbb{E}\[H\_{\pi\pi}\] = -\frac{1}{\pi(1-\pi)}, \qquad
\mathbb{E}\[H\_{\theta\pi}\] = 0, \qquad \mathbb{E}\[H\_{\theta\theta}\]
= (1-\pi)\\\mathbb{E}\[H_W\].\$\$ The hurdle block is the Bernoulli
information, the mixed blocks are exactly zero, and the parent block is
the parent's own expectation weighted by the probability \\1-\pi\\ that
an observation reaches it.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations, passed to the parent, which uses its
  length.

- theta:

  A named list with the parent's parameters followed by `za`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  Ignored for the hurdle and mixed blocks. It is also not forwarded to
  the parent, which takes its own default. Present so that the signature
  matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length `length(y)`, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

`approx` and `nsim` are accepted so that the signature matches the
generic's and neither is read here. They are also NOT forwarded to the
parent, so a parent whose own expected Hessian is approximate takes its
own defaults.

## Notation

\\f_W\\ is the parent's density, \\\pi\\ the probability of the atom at
zero, \\f_Y\\ the mixed density and \\\ell\\ the log-density of one
observation.

## See also

[`distrib_hessian.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedContinuousDistrib.md)
for the observed matrix,
[`distrib_expected_hessian.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedDiscreteDistrib.md)
for the discrete branch, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)
set.seed(4)
y <- distrib_rng(d, 300, theta)

EH <- distrib_expected_hessian(d, y, theta)
vapply(EH, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma       za_za    mu_sigma       mu_za    sigma_za 
#>   -0.175000   -0.350000   -4.761905    0.000000    0.000000    0.000000 

# The hurdle block is the Bernoulli information and the mixed ones vanish.
c(reported = EH$za_za[1], bernoulli = -1 / (0.3 * 0.7),
  mixed = EH$mu_za[1])
#>  reported bernoulli     mixed 
#> -4.761905 -4.761905  0.000000 

# The parent block is the parent's own, weighted by 1 - pi.
c(reported = EH$mu_mu[1], weighted_parent = 0.7 * (-1 / 2^2))
#>        reported weighted_parent 
#>          -0.175          -0.175 
```
