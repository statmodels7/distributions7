# Truncated Random Number Generator (Discrete)

Draws by inverse transform on the parent, exact for a discrete
distribution function: \\Y = Q\\F(L^-) + V Z\\\\ with \\V\\ uniform on
\\(0, 1)\\. Every draw is a support point inside \\\[L, U\]\\, in one
pass and with no rejection.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- n:

  Number of observations to generate, a single positive integer.

- theta:

  A named list of the parent's parameters. A component varying by
  observation must have length `n`.

## Value

A numeric vector of `n` random draws, every one a support point in
\\\[L, U\]\\.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_rng()`](https://statmodels7.github.io/distributions7/reference/trunc_rng.md)
for the shared body,
[`distrib_quantile.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedDiscreteDistrib.md),
which does the work, and
[`distrib_rng.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedContinuousDistrib.md)
for the continuous branch.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

set.seed(4)
y <- distrib_rng(ztp, 5000, theta)
c(zeros = sum(y == 0), min = min(y))
#> zeros   min 
#>     0     1 

# The sample mean tracks the zero-truncated mean, mu / (1 - exp(-mu)).
c(sampled = mean(y), theory = 2 / (1 - exp(-2)))
#>  sampled   theory 
#> 2.304000 2.313035 
```
