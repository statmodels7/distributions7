# Random Starting Values

The base method: `n_start` independent draws from
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md),
which samples each parameter inside its own domain and never looks at
the response. It is what a family gets when it registers no method of
its own, and it is adequate wherever the parameters are of order one
whatever the data, which is true of a shape, a dispersion or a
probability and false of a location or a scale.

## Arguments

- distrib:

  A
  [`distrib()`](https://statmodels7.github.io/distributions7/reference/distrib.md)
  object.

- y:

  The response. Unused here, and accepted only because the generic
  passes it.

- n_start:

  How many starting values to draw, a single positive integer. Defaults
  to 5. A value below 1 is raised to 1.

- ...:

  Unused.

## Value

A list of `n_start` named parameter lists on the parameter scale.

## See also

[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic;
[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md),
which the univariate classes register instead;
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
for the draw itself.
