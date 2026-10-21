# Shared Argument Handling for the Derivative Generics

Aligns `theta`, checks that every parameter has length 1 or \\n\\, and
recycles a scalar `y` up to \\n\\ when `theta` is vectorized. A constant
of the family that varies by observation (a binomial's vector `size`,
including that of its parents) enters \\n\\ too: the response and the
parameters must then have length 1 or the length of that constant, and
any other length signals an error.

## Usage

``` r
check_derivative_args(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A list with elements `y` and `theta`, which are conformable.

## Details

An empty `y` is allowed through untouched, giving empty derivatives, as
`dnorm(numeric(0))` gives `numeric(0)` and as
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
already behaves. Without the special case the recycling check below
rejects it with the nonsensical message
`"'y' must have length 1 or 1, not 0"`.

## See also

[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
[`check_params_dim()`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md)
