# The Mixing Weight Starts at the Proportion of Zeros

Returns the observed proportion of zeros \\\hat p_0\\ as the start of
the zero-inflation probability, kept inside \\\[1/(2n), 1 - 1/(2n)\]\\.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- y:

  The response, a numeric vector.

- ...:

  Unused.

## Value

A named list of one number, named after the mixing parameter.

## Details

Under the model \\P(Y = 0) = \pi + (1 - \pi) f(0) \ge \pi\\, so \\\hat
p_0\\ over-estimates \\\pi\\ and does so from the side away from the
flat end of the chart, which is what makes it a safe start. Measured
against the intercept-only start over three parents (negative binomial,
Poisson-inverse Gaussian, Poisson), with and without inflation, three
samples each: never worse in log-likelihood, better in three of the nine
samples WITHOUT inflation (by 0.43, 0.05 and 3.03, where the default
stopped at the edge and an interior maximum was higher), and converging
in all eighteen. On the sample where the default stopped at \\\pi = 0\\
with an inflated negative binomial, every start with a linear predictor
between -5 and -1.1 reached the interior maximum and every start at -7
or below did not.

The value is on the parameter scale; the layer that uses it carries it
onto the parameter's own link, whichever it is.

## See also

[`distrib_intercept_start()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.md)
for the generic,
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).
