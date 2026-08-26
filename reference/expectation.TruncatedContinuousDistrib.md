# Expectation for Truncated Continuous Distributions

Computes \\\mathbb{E}\[f(Y)\]\\ by integrating the density over \\\[L,
U\]\\ and then adding each surviving atom's contribution. The integral
alone is what the inherited continuous method gives, and it is correct
unless the parent carries point masses, as
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
of a continuous distribution does.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- f:

  A function `f(y, theta, ...)`, vectorized in `y` as the inherited
  method requires.

- theta:

  A named list of the parent's parameters.

- ...:

  Passed to `f` and to the inherited integration.

## Value

A single number, the expectation of `f` under the truncated law.

## Details

The atoms come from
[`distrib_atoms.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md),
already rescaled by \\1/Z\\, so each contributes `p * f(y, theta, ...)`.
With no atoms the loop is empty and the result is exactly the inherited
method's.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
for the generic,
[`distrib_atoms.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md)
for the masses it adds, and
[`expectation.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/expectation.ZeroAdjustedContinuousDistrib.md),
which splits the same way before truncation.

## Examples

``` r
# A mixed parent: the mass at zero is part of the answer.
tz <- truncated(zero_adjusted(gamma2_distrib()), upper = 5)
th <- list(mu = 2, sigma2 = 1, za = 0.3)

# The total mass is one only because the atom is counted.
expectation(tz, function(y, theta) rep(1, length(y)), th)
#> [1] 1

# Which is what mean() and variance() are built on.
c(mean = mean(tz, th), variance = variance(tz, th))
#>     mean variance 
#> 1.368951 1.414969 

# With no atoms it agrees with plain integration over the interval.
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)
fy <- function(y) y * distrib_pdf(tn, y, theta)
c(method = expectation(tn, function(y, theta) y, theta),
  integral = integrate(fy, -1, 2)$value)
#>    method  integral 
#> 0.4159512 0.4159512 
```
