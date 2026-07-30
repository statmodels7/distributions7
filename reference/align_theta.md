# Align Parameters to the Distribution's Parameter Order

Coerces `theta` to a list and, when it is named, reorders it to match
`distrib@params`, so that methods can safely access parameters by
position. This prevents silently wrong results when a user supplies a
named list in a different order (e.g. `list(sigma = 2, mu = 0)` for a
Gaussian).

## Usage

``` r
align_theta(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list (or named numeric vector) of parameters, or an unnamed
  list/vector given in the order of `distrib@params`.

## Value

A list whose first `distrib@n_params` elements correspond to
`distrib@params`, in that order.
