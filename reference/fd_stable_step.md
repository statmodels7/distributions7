# A Parameter Step That Chooses Itself

Evaluates a difference quotient in one parameter at both steps of
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
and returns the one that agrees better with itself at half its own step,
together with the quotient there. Every numerical derivative in the
parameter direction takes its step from it.

## Usage

``` r
fd_stable_step(quotient, theta_j, bounds_j, h_rel, need_value = TRUE)
```

## Arguments

- quotient:

  A function of one step, returning the difference quotient at that
  step: a numeric vector, or a list of them where the caller differences
  several components at once, in which case the reading is taken over
  all of them together.

- theta_j:

  A numeric vector, the values of one parameter.

- bounds_j:

  That parameter's domain, a length-2 numeric vector, or `NULL`.

- h_rel:

  The relative step size.

- need_value:

  Whether the caller will use the quotient at the chosen step. `FALSE`
  for a caller that wants only the step, and then nothing is evaluated
  at all where the two candidates coincide: the higher orders read only
  `h`, and their quotient is a difference of whole Hessians.

## Value

A list of two elements: `h`, the step chosen, and `value`, the quotient
at it, or `NULL` when `need_value` is `FALSE` and no quotient had to be
taken.

## Details

For each of the two steps the quotient \\q(h)\\ is also taken at
\\h/2\\, and the step kept is the one whose
[`fd_self_consistency()`](https://statmodels7.github.io/distributions7/reference/fd_self_consistency.md)
reading is strictly the smaller. Ties, and readings that are not usable,
keep the clamped step, so the answer is the clamped one unless the
scaled step is measurably better.

The choice is made once for the whole vector rather than observation by
observation, which is the opposite of what
[`fd_stable_quotient()`](https://statmodels7.github.io/distributions7/reference/fd_stable_quotient.md)
does in the response direction, and the difference is not an oversight.
There each observation carries its own \\y\\ and therefore its own step,
so a per-observation choice is a choice between two quantities that
genuinely differ; here the parameter is usually one number for the whole
sample, the two candidate quotients estimate the same thing, and the
variation of the self-consistency reading across observations is the
rounding rather than a signal. Measured over a census of 324 cells at
distances from 1 to \\10^{-8}\\ from a bound, the whole-vector choice
puts 307 gradients and 238 diagonal Hessian components within
\\10^{-6}\\ of the analytic value against 285 and 229 for the
per-observation choice, and on the 46 cells where the two differ the
whole-vector one is the better on 44.

It costs four quotients where one would do, except where the two steps
coincide at every entry, which is every call on a parameter with no
finite bound and every positive parameter at or above one: there the
clamped step is returned at once and nothing is evaluated twice.
Measured on the same census, 61 cells of 324 take that path and their
result is [`identical()`](https://rdrr.io/r/base/identical.html) to what
the clamped step alone produced.

Each order chooses its own step, with its own `h_rel` and its own
stencil. A single choice made at first order and reused would cost 2
cells of 324 on the Hessian; the two tests agree on 244 of the 263 cells
where the steps differ at all.

## See also

[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
for the two candidates,
[`fd_self_consistency()`](https://statmodels7.github.io/distributions7/reference/fd_self_consistency.md)
for the reading they are compared on, and
[`fd_stable_quotient()`](https://statmodels7.github.io/distributions7/reference/fd_stable_quotient.md)
for the response direction.
