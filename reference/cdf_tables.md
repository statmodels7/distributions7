# CDF Derivative Tables up to a Given Order

Assembles the derivatives of \\F\\ of orders 1 to `order`. A discrete
family uses the exact finite sum at every order. A continuous family
uses its own
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
at orders 1 and 2 where it registers them, and
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
otherwise; at orders 3 and 4 it uses
[`continuous_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/continuous_cdf_deriv_k.md),
which is the compiled quadrature for some families and one product
stencil on the distribution function for the others.

## Usage

``` r
cdf_tables(distrib, q, theta, order, closed = TRUE)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The highest order wanted, 1 to 4.

- closed:

  Logical of length 1. `TRUE`, the default, takes orders 1 and 2 from
  the family's own
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
  and
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
  where
  [`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
  finds them; `FALSE` takes them from
  [`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md).
  A caller that is itself registered as one of those two methods passes
  `FALSE`, which keeps it from calling itself.

## Value

A list of length `order`. Element \\k\\ is a named list of \\k\\-th
derivatives of \\F\\, on the natural scale and the lower tail, keyed as
[`deriv_names(distrib@params, k)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

Every order below the one wanted is collected, and not just the one
wanted, because the moment-to-cumulant relation that
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md)
applies is a sum over partitions of the multi-index: a partition into
\\k\\ blocks reads \\k\\ lower-order ratios.

Taking orders 1 and 2 from the family's own methods matters on the log
scale, where
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md)
combines them with the higher orders: a differenced table at order 1 or
2 would limit the accuracy of an exact one at order 3 or 4.

## See also

[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md),
the consumer;
[`discrete_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md)
and
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md),
the two routes.
