# CDF Derivatives of a Discrete Distribution at Any Order

The general form of
[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md):
evaluates \\\partial^I F(q) = \sum\_{y \le q} f(y)\\\partial^I
f/f\\(y)\\ for a discrete family at any order up to four. Nothing is
differenced; the sum is exact wherever the support has a finite lower
bound, which the discrete class requires.

## Usage

``` r
discrete_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib`.

- q:

  A numeric vector of quantiles. Each is handled separately, the support
  being walked up to it, so the cost grows with the largest one.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The derivative order, 1 to 4.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

What changes above second order is only the summand. \\\partial^I f/f\\
is the complete Bell polynomial in the derivatives of \\\log f\\, and
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
already computes it for the distribution wrappers, so this function
reuses it and carries no second copy of the enumeration.

## Notation

\\f\\ is the mass function, \\F\\ the distribution function, \\\ell =
\log f\\ and \\\partial^I\\ a derivative with respect to a multi-index
of parameters.

## See also

[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
for orders 1 and 2;
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md),
the continuous route;
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
for the summand;
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md),
the caller.
