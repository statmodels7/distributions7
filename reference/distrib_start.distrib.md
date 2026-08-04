# Random Starting Values

The default: `n_start` draws from
[`generate_random_theta`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md),
which uses the parameter domains and not the data. A family with a
better idea registers its own method.

## Arguments

- distrib:

  A
  [`distrib`](https://statmodels7.github.io/distributions7/reference/distrib.md)
  object.

- y:

  The response, unused here.

- n_start:

  How many to draw.

- ...:

  Unused.

## Value

A list of named parameter lists.
