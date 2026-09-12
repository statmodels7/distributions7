# The Self-Consistency of a Difference Quotient Under Step Halving

Returns \\\max\|q(h/2) - q(h)\| / \max\|q(h)\|\\, the reading
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md)
compares its two candidate steps on. It is `Inf` where the half-step
quotient is not finite and `NA` where the full-step one is identically
zero, neither of which is a usable reading.

## Usage

``` r
fd_self_consistency(x_half, x_full)
```

## Arguments

- x_half:

  The quotient at half the step, a numeric vector.

- x_full:

  The quotient at the full step, the same length.

## Value

A single number, possibly `Inf` or `NA`.

## See also

[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md).
