# Where a Univariate Family Starts, From the Data

A starting value whose location and spread are the response's, for any
univariate family that declares what its parameters mean.

## Usage

``` r
start_from_moments(distrib, y, n_start = 5L, ...)
```

## Arguments

- distrib:

  A univariate distribution.

- y:

  The response.

- n_start:

  How many starting values.

- ...:

  Unused.

## Value

A list of named lists, one per start.

## Details

The base method draws each parameter from its own domain and never looks
at `y`, which is fine while the response is of order one and fails
completely when it is not: on a response of mean 919 and standard
deviation 169 the draws are of order one, the first Newton step is taken
from a point where the residuals are hundreds of standard deviations
wide, and the scale runs to the largest representable double. Measured
on a gaussian,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
recovers N(5, 2) and N(50, 20) and fails on N(500, 200): the defect is a
threshold in the scale of the data, not in the family.

What makes a general fix possible is that every shipped family already
declares `params_interpretation`. A parameter that means a location is
started at the sample median, one that means a spread at the sample
standard deviation or its square, and one whose meaning is a shape, a
dispersion or a probability is left to the draw, those being of order
one whatever the data. A family declaring nothing recognizable loses
nothing: it keeps the draw it had.

The result is CLAMPED strictly inside each parameter's bounds, because a
sample median can land exactly on the boundary of a support and the
validator treats bounds as open. Only the first start is data-based; the
rest stay random, so a caller asking for several still explores.

## See also

[`distrib_start`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
