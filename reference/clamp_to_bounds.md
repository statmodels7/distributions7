# Move a Starting Value Strictly Inside the Parameter Bounds

Returns `theta` with every component moved strictly inside its own
bounds. A moment estimate can land exactly on a boundary. A sample of
non-negative counts gives a median of zero and a sample of proportions a
maximum of one, and
[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md)
treats `params_bounds` as **open**, so such a value is rejected by the
first generic that sees it.

A value at a finite non-zero bound is moved in by 16 machine epsilons
relative to the bound. A value at a bound of exactly zero cannot be,
there being no relative scale there, so it is replaced by \\10^{-8}\\
times the larger of the value's own magnitude and 1. Both are far below
anything a first optimizer step resolves.

## Usage

``` r
clamp_to_bounds(theta, distrib)
```

## Arguments

- theta:

  A named list of parameters on the parameter scale. Components with no
  entry in `params_bounds`, and bounds that are not a pair, are left
  alone.

- distrib:

  The distribution whose `params_bounds` apply.

## Value

`theta`, each component strictly inside its own bounds.

## See also

[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md)
and
[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md),
the callers;
[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
which enforces the open bounds this exists to satisfy.
