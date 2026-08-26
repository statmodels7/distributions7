# Distribution Function of a Truncated Distribution

Evaluates \\F_T(q) = (F(q) - F(L^-))/Z\\, clamped to \\\[0, 1\]\\. One
of the shared method bodies, registered on both truncated classes.

## Usage

``` r
trunc_cdf(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical, default `TRUE`. When `FALSE` the upper tail \\1 - F_T(q)\\ is
  returned.

- log.p:

  Logical, default `FALSE`. When `TRUE` the probability is returned on
  the log scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, of the length `q` and `theta` recycle
to.

## Details

The clamp is there to make the endpoints exact. It returns `0` at and
below \\L\\ and `1` at and above \\U\\, where the unclamped ratio would
give a small negative number or one plus a rounding.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_cdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md)
and
[`distrib_cdf.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedDiscreteDistrib.md),
the two registrations, and
[`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md)
for the inverse.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

# Exactly 0 and 1 at the endpoints.
distrib_cdf(tn, c(-1, 0.3, 2), theta)
#> [1] 0.0000000 0.4609908 1.0000000

# Against the ratio written out.
cs <- distributions7:::trunc_constants(tn, theta)
(pnorm(0.3, 0.3, 1.2) - cs$Fl) / cs$Z
#> [1] 0.4609908

# The two tails sum to one, on either scale.
distrib_cdf(tn, 0.3, theta) + distrib_cdf(tn, 0.3, theta, lower.tail = FALSE)
#> [1] 1
exp(distrib_cdf(tn, 0.3, theta, log.p = TRUE))
#> [1] 0.4609908
```
