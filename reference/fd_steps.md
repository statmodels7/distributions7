# Finite-Difference Steps That Respect a Parameter's Domain

Builds the step \\h\\ for a central difference in one parameter: scaled
by the parameter's magnitude, then shrunk so that \\\theta \pm h\\ stays
strictly inside the parameter's mathematical domain.

## Usage

``` r
fd_steps(theta_j, bounds_j, h_rel)
```

## Arguments

- theta_j:

  A numeric vector, the values of one parameter.

- bounds_j:

  A length-2 numeric vector giving that parameter's domain, or `NULL`.

- h_rel:

  The relative step size, typically a root of machine epsilon chosen for
  the stencil in use.

## Value

A numeric vector of steps, the same length as `theta_j`.

## Details

The domain clamp is what allows a finite-difference fallback to be
offered at all. Parameter domains here are **open** – a scale parameter
is positive, not non-negative – so a step chosen from the magnitude
alone will step a small \\\sigma\\ straight through zero, and the
log-density comes back `NaN` for reasons that look like a bug in the
density. Clamping to 49\\ distance to each finite boundary keeps both
evaluation points inside.

A parameter already on or outside its boundary cannot be rescued this
way, and is reported rather than differentiated.

## See also

[`numerical_gradient`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md),
[`numerical_hessian`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md)
