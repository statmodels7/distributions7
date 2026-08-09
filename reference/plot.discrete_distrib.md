# Plot Method for Discrete Distributions

Visualizes the Probability Mass Function (PMF) of a discrete
distribution object.

## Arguments

- x:

  An object of class `"discrete_distrib"`.

- theta:

  A named list or vector of parameters matching `x@params`. A component
  of length \\k \> 1\\ draws \\k\\ sets of stems.

- xlim:

  Optional numeric vector of length 2 indicating the x-axis range.

- legend:

  Whether to draw the key when several settings are plotted.

- ...:

  Additional arguments passed to the base
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function. `col`
  and `lty` given here win and are recycled over the settings.

## Value

The distribution, invisibly.

## Details

A parameter given as a vector asks for one setting per element, exactly
as in
[`plot.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md).
Several settings are drawn as several sets of stems, shifted sideways so
that one does not stand in front of another, and separated by color and
line type together.
