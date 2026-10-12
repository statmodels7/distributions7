# Skew Normal Mixed Derivatives

Closed form throughout, written in the inverse Mills ratio \\R\\ at \\t
= \alpha z\\. The location and scale components are the location-scale
identity with \\\ell^{(y)} = (\alpha R - z)/\sigma\\ and \\\ell^{(yy)} =
(\alpha^2 R' - 1)/\sigma^2\\ substituted, and the shape component is
\\(R + t R')/\sigma\\, the response entering \\\alpha\\ only through
\\t\\. The three are assembled from a single evaluation of
[`numericals7::mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.html)
rather than by calling
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
which would evaluate it twice more.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- scale:

  Handled by the generic after dispatch; this method returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list with components `mu`, `sigma` and `alpha`, each a numeric
vector of length `length(y)`. All three are exact.

## See also

[`numericals7::mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.html),
evaluated once for all three;
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the family;
[`distrib_cross_y.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.SkewTDistrib.md),
where the shape components are taken from the Hessian.
