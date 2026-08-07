# Response Derivative of a Zero-Adjusted Distribution

Evaluates a response derivative of the parent away from the atom, and
returns `NaN` at it.

## Usage

``` r
za_y_deriv(distrib, y, theta, fun)
```

## Arguments

- distrib:

  A zero-adjusted distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, including the atom probability.

- fun:

  The parent's response-derivative function, e.g.
  [`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md).

## Value

A numeric vector as long as `y`, `NaN` wherever `y == 0`.

## Details

The log-density jumps at zero – \\\log \pi\\ on one side, \\\log((1-\pi)
f(y))\\ on the other – so no derivative in \\y\\ exists there. The
finite-difference default inherited from
[`continuous_distrib`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
would straddle the jump and return a number for it, which is worse than
an error. Away from zero the \\1-\pi\\ factor is constant in \\y\\, so
the parent's own derivative is exact and nothing needs correcting.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
