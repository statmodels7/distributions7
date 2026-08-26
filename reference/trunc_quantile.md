# Quantile Function of a Truncated Distribution

Inverts
[`trunc_cdf()`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
through the PARENT's quantile function, no root-finding of its own being
needed: \\F_T(q) = p\\ holds exactly when \\F(q) = F(L^-) + pZ\\. The
generalized inverse of a discrete distribution function satisfies the
same relation, so the discrete case needs no separate treatment and the
body is registered on both classes.

## Usage

``` r
trunc_quantile(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- p:

  A numeric vector of probabilities, clamped to \\\[0, 1\]\\.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical, default `TRUE`. When `FALSE`, `p` is read as an upper-tail
  probability.

- log.p:

  Logical, default `FALSE`. When `TRUE`, `p` is given on the log scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of the length `p` and `theta` recycle to.

## Details

`p` and `theta` are expanded to a common length before being combined.
\\F(L^-)\\ and \\Z\\ follow the length of `theta` and `p` follows its
own, so multiplying them directly would recycle silently where the two
disagree.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_quantile.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedContinuousDistrib.md)
and
[`distrib_quantile.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedDiscreteDistrib.md),
the two registrations,
[`trunc_cdf()`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
for the inverse, and
[`trunc_rng()`](https://statmodels7.github.io/distributions7/reference/trunc_rng.md),
which draws through it.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_quantile(tn, c(0.1, 0.5, 0.9), theta)
#> [1] -0.6365189  0.3918926  1.5105956

# It inverts the distribution function exactly.
p <- c(0.1, 0.5, 0.9)
max(abs(distrib_cdf(tn, distrib_quantile(tn, p, theta), theta) - p))
#> [1] 2.220446e-16

# And every quantile lies inside the interval.
range(distrib_quantile(tn, seq(0, 1, by = 0.25), theta))
#> [1] -1  2

# A discrete parent needs no separate treatment.
ztp <- truncated(poisson_distrib(), lower = 1)
distrib_quantile(ztp, c(0.1, 0.5, 0.9), list(mu = 2))
#> [1] 1 2 4
```
