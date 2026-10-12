# Which Observations the Finite-Difference Reference Can Be Trusted At

Flags the observations where a central difference has actually
converged, so that
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
compares an analytical derivative only against a reference that is
itself reliable.

## Usage

``` r
fd_is_reliable(fd_at, ref, h_rel, smooth_all)
```

## Arguments

- fd_at:

  A function of one argument, the relative step, returning the
  finite-difference reference computed with that step.

- ref:

  The reference already computed at `h_rel`, a named list of component
  vectors.

- h_rel:

  The relative step `ref` was computed with.

- smooth_all:

  Logical; `TRUE` when every parameter is declared smooth, in which case
  no observation is ever dropped and the guard does not run at all.

## Value

A logical vector as long as the components of `ref`.

## Details

A log-likelihood with a kink has no derivative exactly at the kink, the
Laplace's location being the example the package ships, and a central
difference straddling it returns a number that is simply wrong. An
observation landing within a step of that point therefore makes the
*reference* invalid, not the analytical value being tested, and
comparing against it reports a failure for code that is right. Because
the draws are random, such a failure cannot be predicted.

Rather than hard-code where a kink is, the reference is recomputed with
the step halved: where the two disagree, finite differencing has not
converged and that observation is dropped. For a smooth log-likelihood
nothing is ever dropped, so the check keeps its full strength: a
gradient made 5\\ still caught.

The two estimates are compared **relative to their own magnitude**
rather than against a denominator floored at one: near a kink both are
tiny yet differ by a factor of two, which a floor of one would hide. If
no observation survives, all are kept, because a systematic disagreement
is a real failure and is reported; the guard is meant for a local
disagreement only.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
