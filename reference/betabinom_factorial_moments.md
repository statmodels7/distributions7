# Falling Factorial Moments of the Beta-Binomial

\\E\[Y^{(k)}\] = n^{(k)}\prod\_{j=0}^{k-1}(a+j)/(a+b+j)\\ for \\k = 1,
\dots, 4\\, from which the raw and then the central moments follow.

## Usage

``` r
betabinom_factorial_moments(a, b, n)
```

## Arguments

- a, b:

  The beta shapes.

- n:

  The number of trials.

## Value

A list of the four factorial moments.
