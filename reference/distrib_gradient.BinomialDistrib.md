# Binomial Score

Computes the first derivative of the binomial log-mass with respect to
the probability, one value per observation, in closed form:
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - \dfrac{n -
y}{1-\mu} = \dfrac{y - n\mu}{\mu(1-\mu)}.\$\$ The residual is measured
against \\n\mu\\, the expected count, and divided by \\\mu(1-\mu)\\. Its
sum vanishes at \\\hat\mu = \sum y_i / \sum n_i\\, the pooled
proportion.

On the **link** scale with the default logit the generic's chain rule
gives \\\partial\ell/\partial\eta = y - n\mu\\, the raw residual: the
logit is the canonical link of this family.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- y:

  A numeric vector of counts of successes.

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

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu`, of length
`max(length(y), length(mu), length(distrib@size))`.

## Notation

\\\ell\\ is the log-mass of one observation, \\n\\ the number of trials
and \\\mu \in (0,1)\\ the success probability, with mean \\n\mu\\ and
variance \\n\mu(1-\mu)\\. \\\eta = \log(\mu/(1-\mu))\\ is the log odds.

## See also

[`distrib_hessian.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BinomialDistrib.md)
for the second derivative,
[`distrib_expected_hessian.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BinomialDistrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)
y <- c(0, 4, 10)
th <- list(mu = 0.3)

# The closed form, written out.
all.equal(distrib_gradient(d, y, th)$mu, (y - 10 * 0.3) / (0.3 * 0.7))
#> [1] TRUE

# On the canonical logit link the score is the raw residual.
distrib_gradient(d, y, th, scale = "link")$mu
#> [1] -3  1  7
y - 10 * 0.3
#> [1] -3  1  7

# The summed score vanishes at the pooled proportion.
set.seed(4)
z <- distrib_rng(d, 500, list(mu = 0.3))
sum(distrib_gradient(d, z, list(mu = sum(z) / (500 * 10)))$mu)
#> [1] 4.385381e-14
```
