# CDF Derivative Tables of Every Order Up To One

Assembles the derivatives of \\F\\ of orders 1 to `order` by whichever
route the class uses: the exact finite sum for a discrete family, one
product stencil on the analytic distribution function for a continuous
one. Keeping the choice of route in a single statement is the point of
the function.

## Usage

``` r
cdf_tables(distrib, q, theta, order)
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

## Value

A list of length `order`. Element \\k\\ is a named list of \\k\\-th
derivatives of \\F\\, on the natural scale and the lower tail, keyed as
[`deriv_names(distrib@params, k)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

Every order below the one wanted is collected, and not just the one
wanted, because the moment-to-cumulant relation
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md)
applies is a sum over partitions of the multi-index: a partition into
\\k\\ blocks reads \\k\\ lower-order ratios.

At orders 1 and 2 the tables come from
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
so a family's own closed forms are used where it has them; only orders 3
and 4 reach the routes named above.

## See also

[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md),
the consumer;
[`discrete_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md)
and
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md),
the two routes.
