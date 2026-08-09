# Plot Method for Continuous Distributions

Visualizes the Probability Density Function (PDF) of a continuous
distribution object.

## Arguments

- x:

  An object of class `"continuous_distrib"`.

- theta:

  A named list or vector of parameters matching `x@params`. A component
  of length \\k \> 1\\ draws \\k\\ curves.

- xlim:

  Optional numeric vector of length 2 indicating the x-axis range.

- legend:

  Whether to draw the key when several settings are plotted.

- ...:

  Additional arguments passed to the base
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function. `col`
  and `lty` given here win and are recycled over the curves.

## Value

The distribution, invisibly.

## Details

A parameter given as a vector asks for one curve per element, so
`plot(d, list(mu = 0, sigma = c(1, 2, 4)))` draws three densities that
share a mean and differ in scale. The curves are separated by color and
by line type together, and the parameters that vary are named in a
legend while those held fixed are stated in the title. See
[`plot_settings`](https://statmodels7.github.io/distributions7/reference/plot_settings.md)
for the rule on lengths.

The horizontal range covers every setting: it runs from the smallest
0.5\\ and clamped to the support.
