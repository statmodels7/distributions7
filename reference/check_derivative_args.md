# Shared Argument Handling for the Derivative Generics

Aligns `theta`, checks that every parameter has length 1 or \\n\\, and
recycles a scalar `y` up to \\n\\ when `theta` is vectorised.

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

A list with elements `y` and `theta`, both conformable.

## Details

An empty `y` is allowed through untouched, giving empty derivatives, the
way `dnorm(numeric(0))` gives `numeric(0)` and the way
[`distrib_pdf`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
already behaves. Without the special case the recycling check below
rejects it with the nonsensical message
`"'y' must have length 1 or 1, not 0"`.

## See also

[`align_theta`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
[`check_params_dim`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md)
