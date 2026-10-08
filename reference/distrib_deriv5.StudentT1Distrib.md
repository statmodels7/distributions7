# Student t Fifth-Order Derivatives

Computes the twenty-one distinct fifth derivatives of the location-scale
Student t log-density in \\\mu\\, \\\sigma\\ and \\\nu\\, closed forms
from a compiled kernel written like the lower orders (see
[`distrib_deriv3.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentT1Distrib.md)).

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- scale:

  One of `"parameter"` (the default) or `"link"`; read by the generic.

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of twenty-one numeric vectors named for the multi-index
they carry.

## See also

[`distrib_deriv4.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.StudentT1Distrib.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
length(distrib_deriv5(d, c(-2.5, 0.3), list(mu = 0.4, sigma = 1.2, nu = 5)))
#> [1] 21
```
