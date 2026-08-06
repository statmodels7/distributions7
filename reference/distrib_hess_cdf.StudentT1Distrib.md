# Student t Log-CDF Hessian

Closed form in the location and scale block; the components involving
the degrees of freedom are differenced, having no elementary form.

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

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
