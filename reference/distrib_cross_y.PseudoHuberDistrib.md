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

A named list with one numeric vector per parameter.
