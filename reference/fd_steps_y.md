# Finite-Difference Steps That Respect the Support

The step for a central difference in the response: scaled by \\\|y\|\\
and shrunk so that \\y \pm h\\ stays strictly inside the distribution's
support.

## Usage

``` r
fd_steps_y(y, bounds, h_rel)
```

## Arguments

- y:

  A numeric vector of observations.

- bounds:

  A length-2 numeric vector, the distribution's support.

- h_rel:

  The relative step size.

## Value

A numeric vector of steps, the same length as `y`.

## Details

The response counterpart of
[`fd_steps`](https://statmodels7.github.io/distributions7/reference/fd_steps.md),
and clamped for the same reason: a Gamma observation close to zero would
otherwise be differentiated using a point outside the support, where the
density is not defined.

## See also

[`fd_steps`](https://statmodels7.github.io/distributions7/reference/fd_steps.md),
[`numerical_grad_y`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md)
