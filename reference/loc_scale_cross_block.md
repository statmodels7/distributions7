# The Location and Scale Components of a Mixed Derivative

Returns the two mixed components a location-scale family has in closed
form, from the response derivatives the family already provides. Where
the response enters the log-density only through \\z = (y -
\mu)/\sigma\\, differentiating \\\ell^{(y)} = g'(z)/\sigma\\ gives
\$\$\frac{\partial^2\ell}{\partial y\\\partial\mu} = -\ell^{(yy)},
\qquad \frac{\partial^2\ell}{\partial y\\\partial\sigma} =
-z\\\ell^{(yy)} - \frac{\ell^{(y)}}{\sigma},\$\$ so no new algebra is
needed for either.

## Usage

``` r
loc_scale_cross_block(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from `distrib`, whose first parameter is a
  location and whose second is a scale.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters with the location first and the scale
  second. The order is positional and is not checked.

## Value

An **unnamed** list of two numeric vectors, each of length `length(y)`:
the location component then the scale one. The caller names them.

## Details

The scale formula uses only that \\\sigma\\ is a scale, so it holds when
there is no location at all and \\z = y/\sigma\\. That covers the
exponential's mean, the Weibull's scale and the generalized Pareto's. A
**shape** parameter is not covered: it is written out per family, or
differenced where it has no elementary form.

Both components are recycled to the length of `y` by the `+ 0 * y` in
the body, so a family whose second response derivative is constant still
returns a vector.

## See also

[`loc_scale_cross_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_y.md)
and
[`partial_loc_scale_cross_y()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_cross_y.md),
the two bodies built on this;
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
the two quantities it reads;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
