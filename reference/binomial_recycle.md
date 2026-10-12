# Recycle the Arguments of the Binomial Kernels

Recycles the response, the success probability and the number of trials
to their common length before they reach the compiled kernels, which
read each argument either as a scalar or at the observation's index.
Each argument must have length 1 or the common length.

## Usage

``` r
binomial_recycle(y, mu, size)
```

## Arguments

- y:

  The response, a numeric vector.

- mu:

  The success probability, a numeric vector.

- size:

  The number of trials, a numeric vector.

## Value

A list with components `y`, `mu` and `size`, each of the common length.
