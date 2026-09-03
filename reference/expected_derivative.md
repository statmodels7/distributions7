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
  approx = c("opg", "bartlett", "integrate", "mc"),
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

`"opg"` and `"bartlett"` are two readings of the second Bartlett
identity and are NOT the same computation. The identity equates
\\\mathbb{E}\[\ell^{(ij)}\]\\ with
\\-\mathbb{E}\[\ell^{(i)}\ell^{(j)}\]\\; `"bartlett"` evaluates the
expectation, which is a sum over the support or a quadrature, while
`"opg"` reads the integrand at the observation and takes no expectation
at all. They agree in expectation and differ by orders of magnitude in
cost, which is why `"opg"` is the default at order 2.

Above order 2 the outer product of scores is not an identity for
anything, so `"opg"` there is routed to `"bartlett"` rather than
refused: the argument is one setting for a whole call and an order-4
method may ask for the same strategy an order-2 one was given.

## See also

[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
