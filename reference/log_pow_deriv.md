# Derivatives of log(p) and log(1 - p)

The two elementary logarithmic derivatives the zero wrappers need: \\d^k
\log p = (-1)^{k-1}(k-1)!/p^k\\ and \\d^k \log(1-p) = -(k-1)!/(1-p)^k\\.

## Usage

``` r
log_pow_deriv(p, k, complement = FALSE)
```

## Arguments

- p:

  A numeric vector of probabilities.

- k:

  The derivative order.

- complement:

  Logical; `TRUE` for \\\log(1-p)\\.

## Value

A numeric vector.
