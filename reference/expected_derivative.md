# Dispatch an Expected-Derivative Strategy

Shared entry point for the fallback expected-derivative methods:
validates `approx` and routes to the chosen strategy.

## Usage

``` r
expected_derivative(
  distrib,
  y,
  theta,
  order,
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000
)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

- approx:

  One of `"bartlett"`, `"integrate"`, `"mc"` or `"opg"`.

- nsim:

  Number of draws, used only by `"mc"`.

## Value

A named list of expected derivative component vectors.

## Details

`"opg"` is folded into `"bartlett"` here rather than implemented
separately, because at order 2 the outer product of gradients *is* the
Bartlett identity; the two names describe the same computation from
different traditions.

## See also

[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
