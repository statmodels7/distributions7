# Student t Log-CDF Gradient

Closed form in the location and scale, \\-f(y)\\ and \\-z f(y)\\; the
degrees of freedom are differenced, having no elementary form.

## Arguments

- distrib:

  A `StudentT1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma` and `nu`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
