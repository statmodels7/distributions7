# Numerical Summation of Discrete Series

Calculates the sum of a function `f(x)` over a sequence of integers from
`start` to `end`. The function is designed to handle finite sums,
one-sided infinite series, and doubly infinite series by automatically
adapting its summation strategy.

## Usage

``` r
numerical_series(
  f,
  start = 0,
  end = Inf,
  step = 10000,
  tol = 1e-10,
  maxit = 1000000L,
  reltol = TRUE
)
```

## Arguments

- f:

  A function taking a vector of integers `x` and returning a vector of
  numeric values. **Must be vectorized**.

- start:

  Numeric. Starting value. Can be finite, `Inf`, or `-Inf`. Defaults to
  `0`.

- end:

  Numeric. Ending value. Can be finite, `Inf`, or `-Inf`. Defaults to
  `Inf`.

- step:

  Integer. Number of terms to calculate in a single vectorized batch.
  Defaults to `1000`.

- tol:

  Numeric. Tolerance threshold for convergence. Defaults to `1e-10`.

- maxit:

  Integer. Safety limit for the maximum number of batch iterations.
  Defaults to `1000000`.

- reltol:

  Logical. If `TRUE` (default), uses a hybrid relative tolerance.

## Value

A numeric scalar representing the calculated sum.

## Details

**1. Summation Strategies:** The function automatically detects the
domain topology based on `start` and `end`:

- **Forward (Standard):** If `start <= end` (e.g., `1` to `Inf`).

- **Backward (Reflection):** If `start > end` (e.g., `-1` to `-Inf`),
  evaluates `f(-x)`.

- **Doubly Infinite (Folding):** If `start == -Inf` and `end == Inf`,
  folds around 0.

**2. Speed and Convergence:** It monitors convergence using the sum of
absolute values in the current chunk, preventing premature stops on
alternating series while maintaining high precision.

**3. Underflow & Divergence Detection:** Includes heuristics to stop
early if the sequence starts growing in absolute terms (divergence), or
skips up to 50 empty chunks to protect from premature stopping when
`f(x)` evaluates exactly to `0` at the start.

## Examples

``` r
# the mean of a Poisson, summed over its support
numerical_series(function(y) y * distrib_pdf(poisson_distrib(), y, list(mu = 2)))
#> [1] 2
```
