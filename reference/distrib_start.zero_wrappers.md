# Starting Values for the Zero Wrappers, Read Off the Data

Starting values for
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).
The mixing probability starts at the observed proportion of zeros \\\hat
p_0\\, kept inside \\\[1/(2n), 1 - 1/(2n)\]\\, and the parent's
parameters start where
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
of the PARENT puts them: on the non-zero observations for a
zero-adjusted family, whose positive part is the parent away from zero,
and on every observation for a zero-inflated one, whose parent also
produces zeros.

## Usage

``` r
start_zero_wrapper(distrib, y, n_start = 5L, keep)
```

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib`, `ZeroAdjustedContinuousDistrib` or
  `ZeroInflatedDistrib`.

- y:

  The response, a numeric vector.

- n_start:

  How many starting values, a single positive integer.

- keep:

  A function of the response returning which observations the parent's
  starts are read from.

## Value

A list of `n_start` named parameter lists on the parameter scale, in the
order of `distrib@params`.

## Details

Without these methods the wrappers took the univariate fallback, which
reads the location off the median of the whole sample: on data with half
their values at zero that is zero, and a start at a mean of zero is the
edge of a log link. On a zero-adjusted Poisson-inverse Gaussian sample
the first start then converged to a degenerate point (\\\mu \to 0\\,
\\\sigma \to \infty\\) 97 log-likelihood units below the maximum.

For a zero-adjusted family \\\hat p_0\\ is the maximum likelihood
estimate of the mixing probability, the likelihood factorizing into a
binomial part and a positive part. For a zero-inflated one it
over-estimates \\\pi\\, since \\P(Y = 0) = \pi + (1-\pi) f(0) \ge \pi\\,
and does so from the side away from the flat end of the chart; see
[`distrib_intercept_start.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.ZeroInflatedDistrib.md).
Only the first start carries \\\hat p_0\\; the others keep a random
mixing probability, and the parent's own starts in the same order.

## See also

[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic,
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
[`distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.md).
