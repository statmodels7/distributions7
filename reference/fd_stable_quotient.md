# A Response Difference That Chooses Its Step

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
  difference quotient at those steps. It should divide by the steps its
  evaluation points actually lie at, as
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

Neither step is right everywhere, which is why both are read. Where the
log-density is singular at a bound the clamped step has an error that
stops falling once the clamp binds, and where it is smooth at the bound
the scaled step is the worse one at the smallest distances. Measured
over
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)'s
response row for 32 continuous families, three parameter values and five
seeds, the clamped reference failed 35 rows of 480, all of families
singular at a bound, and this one fails none, its worst at 5.0e-05. On a
sweep to within \\10^{-8}\\ of a bound it has no failure among 264
points of the singular families, where the clamped step fails 165 on the
first derivative and 198 on the second, and it fails 1 and 47 of 144 on
the families smooth at the bound, where the clamped step fails 0 and 43.
On the distribution function the same choice took the check's grid from
3 failed rows to none.

Where it keeps the scaled step and the clamped one would have passed,
which is 6 points of that sweep's 816, the scaled numerator is a whole
number of units in the last place of the log-density that scales exactly
with the step, so the quotient agrees with itself to the bit at half the
step and the test reads it as perfectly stable.

It costs four quotients where one would do, except where the two steps
coincide at every observation, which is every call on a family with no
finite bound: there the clamped quotient is returned at once.

## Near a bound that is not zero

Near zero the spacing of doubles is relative, so a step scaled on the
distance can be as small as the distance asks. Near any other bound the
spacing is absolute, one unit of the last place at the bound, and within
about \\10^{-10}\\ of it two things go wrong. A step that is not a whole
number of units is not the step the evaluation points lie at, so a
quotient dividing by the nominal step is out by their rounding; and the
scaled step falls below one unit and rounds to zero. The choice cannot
see either, reading the same arithmetic at both steps. So no candidate
step is shorter than \\2\lvert y\rvert\varepsilon\\, which leaves its
half step at least one unit long, and the quotient a caller passes
divides by the steps actually taken, as
[`fd_first_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
and
[`fd_second_taken()`](https://statmodels7.github.io/distributions7/reference/fd_taken.md)
do. Measured on
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
at `mu = 0.5, phi = 0.5`, singular at 1, the response row of 100 draws
failed in 486 samples of 2000 with the nominal quotient and fails in 10
with both repairs; over the census above no verdict moves and the worst
statistic goes from 3.0e-05 to 5.0e-05, and on the distribution
function's grid no statistic moves by more than 1.4e-10. What neither
repair reaches is a point within a few units of such a bound, where the
relative error of the reference falls as about \\5.4/d^2\\ in the
distance \\d\\ counted in units;
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
