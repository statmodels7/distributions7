# Expand Parameters to Common Length

Expands scalar parameters in a list to match the maximum length found
(or a specified length), ensuring all vectors are ready for element-wise
operations.

## Usage

``` r
expand_params(theta, n)
```

## Arguments

- theta:

  A named list of parameters.

- n:

  (Optional) The target length. If missing, defaults to
  `max(lengths(theta))`.

## Value

A list where all elements have length `n`.
