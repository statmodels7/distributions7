# A Response Difference with a Selected Step

Evaluates a central difference in the response at both steps of
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
and keeps, observation by observation, the one that agrees better with
itself at half its step. It is the reference
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
compares a family's response derivatives and its distribution function
against.

## Usage

``` r
fd_stable_quotient(quotient, y, bounds, h_rel)
```

## Arguments

- quotient:

  A function of a vector of steps as long as `y`, returning the
  difference quotient at those steps. It should divide by the steps at
  which its evaluation points actually lie, as
  [`fd_first_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
  does.

- y:

  The points, a numeric vector.

- bounds:

  The support, a numeric vector of length two.

- h_rel:

  The relative step, a single positive number.

## Value

A numeric vector as long as `y`.

## Details

For each of the two steps the quotient \\q(h)\\ is also taken at
\\h/2\\, and its self-consistency is \\\|q(h/2) - q(h)\|/\|q(h)\|\\:
infinite where \\q(h/2)\\ is not finite, and undefined where \\q(h)\\ is
exactly zero. The step scaled on the distance to the bound is kept only
where both readings are defined and its own is strictly the smaller;
everywhere else, ties included, the clamped step is kept. The value
returned is the quotient at the full step.

Neither step is accurate for every family, which is why both are read.
Where the log-density is singular at a bound, the error of the clamped
step stops falling once the clamp binds; where it is smooth at the
bound, the scaled step is the less accurate one at the smallest
distances. Where both readings are defined, the clamped step is kept
unless the scaled one agrees strictly better with itself at half the
step.

In a few cases the scaled step is kept although the clamped step would
have been accurate. The scaled numerator is then a whole number of units
in the last place of the log-density, which scales exactly with the
step, so the quotient agrees with itself at half the step to the last
bit and the test reads it as stable.

It costs four quotients where one would do, except where the two steps
coincide at every observation, which is every call on a family with no
finite bound: there the clamped quotient is returned at once.

## Near a bound that is not zero

Near zero the spacing of doubles is relative, so a step scaled on the
distance can be as small as the distance requires. Near any other bound
the spacing is absolute, one unit of the last place at the bound, and
within about \\10^{-10}\\ of it two things go wrong. A step that is not
a whole number of units is not the step at which the evaluation points
lie, so a quotient that divides by the nominal step is out by their
rounding; and the scaled step falls below one unit and rounds to zero.
The selection rule cannot tell these cases apart, because both
candidates come from the same arithmetic. Every candidate step is
therefore at least \\2\lvert y\rvert\varepsilon\\, which leaves its half
step at least one unit long, and the quotient that a caller passes
divides by the steps actually taken, as
[`fd_first_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
and
[`fd_second_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
do. Neither repair reaches a point within a few units of such a bound,
where the relative error of the reference decreases as about \\5.4/d^2\\
in the distance \\d\\ counted in units;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
leaves those draws out.

## See also

[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
for the two steps,
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
for where it is used,
[`fd_first_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
and
[`fd_second_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
for quotients on the steps taken.
