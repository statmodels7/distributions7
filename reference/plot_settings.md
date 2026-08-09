# The Parameter Settings a Plot Draws

Splits a `theta` whose components may be vectors into one setting per
curve, and reports which parameters vary across them.

## Usage

``` r
plot_settings(x, theta)
```

## Arguments

- x:

  A distribution object.

- theta:

  A named list or numeric vector of parameters.

## Value

A list with `settings` (one scalar `theta` per curve), `k` and `varying`
(the names of the parameters that differ between settings).

## Details

A component of length \\k \> 1\\ asks for \\k\\ curves and a component
of length one is held across all of them, so
`list(mu = 0, sigma = c(1, 2, 4))` is three settings sharing a mean.
Every component must therefore have length one or the same \\k\\;
anything else is rejected rather than recycled, since a length that
divides \\k\\ is far more likely to be a mistake than a request.

This meaning is available because a plot has no data to recycle against.
The derivative and density generics read a vector component as one value
per observation, which is a different question asked of the same object.
