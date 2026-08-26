# Mixed Derivatives of a Reparametrized Distribution

The parent's mixed block carried by the first-order chain rule on the
map.

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of the new parameters.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter of the **new**
parametrization, keyed by `distrib@params`, each of length `length(y)`.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
for the wrapper;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic;
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
where the chain rule is the same shape.
