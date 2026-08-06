# Generalized Pareto Log-CDF Gradient

Closed form. With \\t = 1 + \xi q/\sigma\\ and \\S = t^{-1/\xi}\\,
\\\partial F/\partial\sigma = -S q/(\sigma^2 t)\\ and \\\partial
F/\partial\xi = -S(\log t/\xi^2 - q/(\xi\sigma t))\\.

## Arguments

- distrib:

  A `GPDDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `sigma` and `xi`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## Details

The shape direction is a difference of two quantities of size \\z/\xi\\
that cancel as \\\xi \to 0\\, so below \\10^{-4}\\ it is taken from the
series \\z^2/2 - 2\xi z^3/3 + 3\xi^2 z^4/4\\ that the cancellation
leaves, matching the branch the distribution function itself takes
there.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
