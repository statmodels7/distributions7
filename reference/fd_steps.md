# Finite-Difference Steps That Respect a Parameter's Domain

Builds the step \\h\\ for a central difference in one parameter: scaled
by the parameter's magnitude and kept strictly inside the parameter's
mathematical domain, either by cutting it to under half the distance to
the nearest finite bound or by scaling it on that distance.

## Usage

``` r
fd_steps(
  theta_j,
  bounds_j,
  h_rel,
  to_bound = c("clamp", "scale"),
  order = 1L,
  accuracy = 2L
)
```

## Arguments

- theta_j:

  A numeric vector, the values of one parameter.

- bounds_j:

  A length-2 numeric vector giving that parameter's domain, or `NULL`.

- h_rel:

  The relative step size, typically a root of machine epsilon chosen for
  the stencil in use.

- to_bound:

  `"clamp"`, the default, or `"scale"`, as above.

- order, accuracy:

  The stencil the step is for, which fix its reach through
  [`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html).
  The defaults are the central difference, whose reach is one.

## Value

A numeric vector of steps, the same length as `theta_j`.

## Details

The domain clamp is what allows a finite-difference fallback to be
offered at all. Parameter domains here are **open**: a scale parameter
is positive, not non-negative. A step chosen from the magnitude alone
therefore takes a small \\\sigma\\ straight through zero, and the
log-density comes back `NaN` for reasons that look like a bug in the
density. Clamping to 49\\ distance to each finite boundary keeps both
evaluation points inside.

`to_bound = "scale"` takes \\h = h\_{rel}\min(\max(1,\|\theta\|), d)\\,
with \\d\\ the distance to the nearest finite bound, so the step shrinks
in proportion as the parameter approaches the bound rather than sitting
at \\0.49\\d\\. Where the log-density is singular in that parameter at
the bound the clamped step's error stops falling once the clamp binds,
and where the log-density instead saturates the scaled step is the worse
of the two, its shorter step letting the rounding dominate. Neither is
right everywhere, and
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md)
is what chooses between them.

The two give the same step wherever \\d \ge \max(1, \|\theta\|)\\, which
is every parameter of a family with no finite bound, and also a positive
parameter at or above one.

A parameter already on or outside its boundary cannot be rescued this
way, and is reported rather than differentiated.

The clamp is on the stencil's OUTERMOST node and not on the step, so it
divides by the reach, the largest offset the rule evaluates at in units
of \\h\\. That number is the stencil's and is read from
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
rather than assumed: a central difference and a three-point second
difference both reach one step out, which is every use this package
makes of the clamped branch today, while the five-point rules of
accuracy four reach two. Without the division a step cut to \\0.49 d\\
puts the outermost node of a reach-two rule at \\0.98 d\\, which is
inside the domain and not what the clamp says it is.

With the reach accounted for, the clamped branch at the default `h_rel`
is
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
with `bounds`, and a test asserts the two are
[`identical()`](https://rdrr.io/r/base/identical.html) so they cannot
drift. It is written here rather than delegated because `h_rel` is an
argument the callers vary:
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md)
halves it and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
reads the quotient at two steps, neither of which `fd_step()` can
express.

The scaled branch needs no such division. Its step is at most \\h\_{rel}
d\\, so the outermost node sits at \\d(1 - r\\h\_{rel})\\ for a reach
\\r\\, inside the domain for any stencil at any root of machine epsilon.

## See also

[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md),
which chooses between the two, and
[`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md)
and
[`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md),
which take the chosen step.
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
is the response counterpart.
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
is the same rule at the default `h_rel`.
