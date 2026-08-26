# Zero-Adjusted Distribution Object

Makes the probability of a zero a parameter of its own, \\\pi\\, carried
by `za`, and leaves everything else to the parent. What that means
depends on the parent's type, and the constructor dispatches on it:

- a DISCRETE parent with 0 in its support gives a HURDLE model. The mass
  the parent puts at zero is removed, the parent is renormalized over
  the positive values, and \\\pi\\ takes its place.

- a CONTINUOUS parent gives a MIXED distribution. Nothing has to be
  removed, \\P(Y = 0)\\ being zero already; a point mass \\\pi\\ is
  placed at zero and the density is scaled by \\1-\pi\\.

It is the wrapper to reach for when zeros come from their own mechanism,
no claim filed or no rainfall, and the parent describes only what
happens once that mechanism has been passed. Where the zeros come partly
from the parent and partly from a separate source, so that no single
zero can be attributed,
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
is the model.

## Usage

``` r
zero_adjusted(distrib, link_za = logit_link())
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` with 0 in its support, or
  from `continuous_distrib`. An already-wrapped parent and a discrete
  support of fewer than `n_params + 2` points are rejected with an error
  saying which condition failed.

- link_za:

  The link carrying \\\pi\\ to the unconstrained scale, a
  [`linkfunctions7::link`](https://statmodels7.github.io/linkfunctions7/reference/link.html)
  object. Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  which keeps it strictly inside \\(0, 1)\\ at every point of the free
  scale.

## Value

An S7 object of class
[ZeroAdjustedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedDiscreteDistrib.md)
or
[ZeroAdjustedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedContinuousDistrib.md),
matching the parent's branch. Its `params` are the parent's followed by
`za`, whose bound is \\(0, 1)\\ and whose interpretation is
`"prob. of zero"`; `distrib_name` is `"zero-adjusted "` followed by the
parent's.

## A discrete parent: the hurdle

\$\$P(Y = y; \theta, \pi) = \begin{cases} \pi & y = 0 \\ (1 -
\pi)\dfrac{f(y; \theta)}{1 - f(0; \theta)} & y \> 0. \end{cases}\$\$ The
division by \\1 - f(0;\theta)\\ is the TRUNCATION, and it is what
separates this from zero-inflation: it makes \\\theta\\ the parameters
of a law on the positive integers instead of the original count process.
The log-likelihood separates into a Bernoulli part in \\\pi\\ and a
truncated part in \\\theta\\, so every mixed block of the information is
exactly zero and the two halves could be fitted separately.

## A continuous parent: the mixed law

\$\$f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi) f_W(y;\theta) \quad (y \ne
0).\$\$ Here \\f_Y\\ is a density against Lebesgue measure plus a point
mass, so it integrates to \\1 - \pi\\; the remainder is the atom, which
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
reports and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
accounts for. The classical members have a parent on \\(0, \infty)\\ or
\\(0, 1)\\, where zero sits on the boundary: the zero-adjusted gamma,
inverse gaussian or lognormal for semicontinuous data, and the
zero-adjusted beta for proportions.

## Choosing between the two wrappers

Zero-inflation can only ADD zeros, \\P(Y = 0) = \zeta + (1-\zeta)f(0)\\
exceeding \\f(0)\\; the hurdle replaces \\f(0)\\ outright and so also
covers FEWER zeros than the parent implies. Where both apply they are
not nested, and they differ in interpretation more than in fit:
zero-inflation keeps \\\theta\\ as the parameters of the original count
process, the hurdle re-reads them as those of a truncated one. Prefer
the hurdle when a zero is observable evidence of a distinct decision,
and zero-inflation when the two kinds of zero are genuinely
indistinguishable.

## What the parent must be

A discrete parent must have 0 in its support: with \\f(0) = 0\\ there is
no mass to remove. Construction also fails where the result would not be
identified:

- the parent already models a probability of zero. Zero-truncating a
  zero-inflated or zero-adjusted parent cancels its zero parameter out
  of the likelihood entirely, leaving an identically zero score;

- the support is too small to carry one more parameter. A distribution
  on \\k\\ points has \\k-1\\ free probabilities, so at least
  `n_params + 2` support points are needed. Zero-adjusting a Bernoulli
  leaves the truncated part concentrated on \\\\1\\\\ and `mu` vanishes
  from the likelihood.

A continuous parent supported on the whole line is accepted and gives a
spike-at-zero model, but \\y = 0\\ then no longer identifies its own
mechanism: the atom sits where the density is positive. A continuous
parent whose support does not reach zero at all, say \\(2, 5)\\, is
accepted with a WARNING, the atom then being disconnected from the rest
of the law, which is legitimate and rarely intended.

## Notation

\\f\\ is a discrete parent's mass function, \\f_W\\ a continuous one's
density, \\\pi\\ the probability of a zero, \\\zeta\\ the inflation
probability of the other wrapper, and \\k\\ the number of support
points.

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
for the mixture counterpart,
[ZeroAdjustedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedDiscreteDistrib.md)
and
[ZeroAdjustedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedContinuousDistrib.md)
for the two classes,
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for the atom the continuous branch declares, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate the result.

## Examples

``` r
# Hurdle Poisson: the mass at zero is exactly za, not dpois(0, mu).
zap <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.3)
distrib_pdf(zap, 0:5, theta)
#> [1] 0.3000000 0.1100310 0.1650464 0.1650464 0.1237848 0.0742709
c(at_zero = distrib_pdf(zap, 0, theta), parent = dpois(0, 3))
#>    at_zero     parent 
#> 0.30000000 0.04978707 

# And it can be BELOW the parent's, which zero-inflation cannot reach.
c(hurdle = distrib_pdf(zap, 0, list(mu = 3, za = 0.01)),
  parent = dpois(0, 3))
#>     hurdle     parent 
#> 0.01000000 0.04978707 

# The likelihood separates, so every mixed block of the Hessian is zero.
set.seed(2)
y <- distrib_rng(zap, 300, theta)
all(distrib_hessian(zap, y, theta)$mu_za == 0)
#> [1] TRUE

# Semicontinuous data: a spike at zero and a gamma above it.
zag <- zero_adjusted(gamma2_distrib())
distrib_atoms(zag, list(mu = 2, sigma2 = 1, za = 0.3))
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 

# A fit recovers both halves, each from its own part of the data.
set.seed(5)
yg <- distrib_rng(zag, 2000, list(mu = 2, sigma2 = 1, za = 0.3))
round(coef(fit_distrib(zag, yg)), 3)
#>     mu sigma2     za 
#>  2.029  1.016  0.302 

# Two refusals, each naming the condition that failed.
try(zero_adjusted(bernoulli_distrib()))
#> Error : zero_adjusted() cannot wrap 'bernoulli': its support has 2 points, so the family has 1 free
#>   probability, while the wrapped distribution would have 2 parameters. They are not
#>   identified -- different parameter values give exactly the same distribution.
#>   A support of at least 3 points is required.
try(zero_adjusted(zero_inflated(poisson_distrib())))
#> Error : zero_adjusted() cannot wrap 'zero-inflated poisson', which already models the probability of a zero.
#>   Stacking the two leaves only their combination identified: the second
#>   parameter has an identically zero score, and any optimizer will wander
#>   along that ridge. Apply exactly one zero wrapper to a plain distribution.
```
