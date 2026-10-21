# Laplace Mixed Derivatives

The mixed derivatives are in closed form from the location-scale
identity, away from the kink at \\y = \mu\\, where the density has no
derivative in the response and the quantity does not exist.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Handled by the generic after dispatch; this method returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list with components `mu` and `sigma`; each is a numeric vector
of length `length(y)`. The `mu` component is exactly zero everywhere,
because the second response derivative of a Laplace is zero almost
everywhere, and the `sigma` component is
\\\mathrm{sign}(y-\mu)/\sigma^2\\, which jumps sign at the kink.

## See also

[`loc_scale_cross_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_y.md),
the shared body;
[`distrib_cross_y.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Laplace2Distrib.md)
for the rate chart;
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
whose value of zero here gives the zero `mu` component.
