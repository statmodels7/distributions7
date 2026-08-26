# Bernoulli Score

Computes the first derivative of the Bernoulli log-mass with respect to
the probability, one value per observation, in closed form:
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} -
\dfrac{1-y}{1-\mu} = \dfrac{y - \mu}{\mu(1-\mu)}.\$\$ The residual is
divided by the variance, so an observation near a probability of 0 or 1
that goes the other way contributes an arbitrarily large score.

On the **link** scale with the default logit the generic's chain rule
gives \\\partial\ell/\partial\eta = y - \mu\\, the raw residual: the
logit is the canonical link of this family.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- y:

  A numeric vector of zeros and ones.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(0, 1)\\.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of one numeric vector, `mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell\\ is the log-mass of one observation and \\\mu \in (0,1)\\ the
success probability, with variance \\\mu(1-\mu)\\. \\\eta =
\log(\mu/(1-\mu))\\ is the parameter on the link scale, the log odds. A
**canonical** link makes the log-mass linear in the sufficient statistic
times \\\eta\\.

## See also

[`distrib_hessian.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BernoulliDistrib.md)
for the second derivative,
[`distrib_expected_hessian.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BernoulliDistrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()
y <- c(0, 1, 1)
th <- list(mu = 0.3)

# The closed form, written out.
all.equal(distrib_gradient(d, y, th)$mu, (y - 0.3) / (0.3 * 0.7))
#> [1] TRUE

# On the canonical logit link the score is the raw residual.
distrib_gradient(d, y, th, scale = "link")$mu
#> [1] -0.3  0.7  0.7
y - 0.3
#> [1] -0.3  0.7  0.7

# The summed score vanishes at the sample proportion, which is the estimate.
set.seed(4)
z <- distrib_rng(d, 1000, list(mu = 0.3))
sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
#> [1] 2.549072e-13
```
