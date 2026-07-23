# Check Consistency of Parameter Dimensions

Validates that all elements in the provided parameter list have
compatible lengths. Each parameter must have a length of either 1
(scalar) or exactly equal to `n`. This ensures safe vector recycling and
dimensional consistency.

## Usage

``` r
check_params_dim(theta, n)
```

## Arguments

- theta:

  A named list of vectors (parameters). Each element represents a
  parameter of a distribution (e.g., `mu`, `sigma`).

- n:

  (Optional) An integer specifying the required maximum length. If not
  provided, it defaults to the maximum length found among the elements
  of `theta`. Providing this argument allows validation against an
  external dimension (e.g., sample size `n`).

## Value

Returns `NULL` invisibly if the check passes.

## Errors

The function throws an error (`stop`) if it detects any parameter with a
length that is neither 1 nor `n`. The error message lists the specific
parameters causing the mismatch.

## Examples

``` r
# --- Case 1: Implicit max length ---
# Valid: all scalars
check_params_dim(list(mu = 1, sigma = 2))

# Valid: mixing scalar and vector
check_params_dim(list(mu = 1:5, sigma = 1))

# Invalid: incompatible lengths (2 vs 3)
if (FALSE) { # \dontrun{
check_params_dim(list(mu = 1:2, sigma = 1:3))
} # }

# --- Case 2: Explicit n ---
# Valid: vector matches n (5)
check_params_dim(list(mu = 1:5, sigma = 1), n = 5)

# Invalid: vector length (3) does not match required n (5)
# This is useful to enforce consistency with a dataset size n = 5
if (FALSE) { # \dontrun{
check_params_dim(list(mu = 1:3, sigma = 1), n = 5)
} # }
```
