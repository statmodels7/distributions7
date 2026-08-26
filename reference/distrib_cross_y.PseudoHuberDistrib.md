# Pseudo-Huber Mixed Derivatives

Closed form throughout. The location and scale components come from the
location-scale identity; the response enters the shape only through \\D
= \sqrt{\nu + r^2/\sigma^2}\\, so the shape component is \\r/(2\sigma^2
D^3)\\ with \\r = y - \mu\\.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `nu`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with components `mu`, `sigma` and `nu`, each a numeric
vector of length `length(y)`. All three are exact; measured against
Richardson on the analytic response gradient the worst is
\\1.3\times10^{-11}\\ relative.

## See also

[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
for the first two components;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for the family;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
