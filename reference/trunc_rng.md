# Random Generation From a Truncated Distribution

Draws by inverse transform through
[`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md),
\\Y = Q(F(L^-) + V Z)\\ with \\V\\ uniform on \\(0, 1)\\. The draw is
exact and terminates in one pass however small \\Z\\ is, where rejection
sampling from the parent would need \\1/Z\\ draws on average. One of the
shared bodies, registered on both classes.

## Usage

``` r
trunc_rng(distrib, n, theta, ...)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- n:

  The number of draws, a single positive integer.

- theta:

  A named list of the parent's parameters. A component varying by
  observation must have length `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws, every one inside \\\[L, U\]\\.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_rng.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedContinuousDistrib.md)
and
[`distrib_rng.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedDiscreteDistrib.md),
the two registrations, and
[`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md),
which does the work.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

set.seed(1)
y <- distrib_rng(tn, 5000, theta)
c(inside = all(y >= -1 & y <= 2), min = min(y), max = max(y))
#>     inside        min        max 
#>  1.0000000 -0.9991523  1.9995546 

# The sample mean tracks the truncated mean, not the parent's.
c(sampled = mean(y), truncated = mean(tn, theta), parent = 0.3)
#>   sampled truncated    parent 
#> 0.4099989 0.4159512 0.3000000 

# An interval carrying little mass costs no more than any other.
far <- truncated(gaussian1_distrib(), lower = 4, upper = 5)
set.seed(2)
range(distrib_rng(far, 1000, theta))
#> [1] 4.000226 4.995428
```
