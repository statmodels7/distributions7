# Expectation for Zero-Adjusted Continuous Distributions

Computes \\\mathbb{E}\[f(Y)\]\\ by splitting the expectation at the
atom, \$\$\mathbb{E}\[f(Y)\] = \pi\\ f(0) + (1-\pi)\\
\mathbb{E}\_W\[f(W)\],\$\$ with the second term the PARENT's own
expectation. Plain numerical integration over the mixed density would
miss the point mass entirely and return \\(1-\pi)\mathbb{E}\_W\[f\]\\,
silently.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- f:

  A function `f(y, theta, ...)`, receiving the full `theta` including
  `za`. It must be vectorized in `y`, as the parent's own
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
  requires.

- theta:

  A named list with the parent's parameters followed by `za`.

- ...:

  Passed to `f` and to the parent's
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Value

A single number, the expectation of `f` under the mixed law.

## Details

`f` receives the FULL `theta`, `za` included, even though the integral
runs over the parent. The method re-attaches the atom probability inside
a wrapper, so a caller's `f` sees the same parameter list at the atom
and away from it.

## Notation

\\\pi\\ is the probability of the atom at zero, \\W\\ the parent's
variable and \\Y\\ the zero-adjusted one.

## See also

[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
for the generic,
[`distrib_atoms.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md)
for the declaration this rests on, and
[`base::mean()`](https://rdrr.io/r/base/mean.html) and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
the two moments built on it.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

# The first two moments, against the split written out.
c(computed = expectation(d, function(y, theta) y, theta),
  theory = 0.7 * 1)
#> computed   theory 
#>      0.7      0.7 
c(computed = expectation(d, function(y, theta) y^2, theta),
  theory = 0.7 * (1^2 + 2^2))
#> computed   theory 
#>      3.5      3.5 

# mean() and variance() are built on the same split.
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>     0.70     3.01 
c(theory_var = 0.7 * 5 - (0.7 * 1)^2)
#> theory_var 
#>       3.01 

# The parent's own expectation misses the atom, and is the wrong answer for
# the mixed law by exactly the factor 1 - pi.
expectation(gaussian1_distrib(), function(y, theta) y,
            theta[c("mu", "sigma")])
#> [1] 1
```
