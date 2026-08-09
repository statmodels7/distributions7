# CDF Derivative Tables of Every Order Up To One

The derivatives of \\F\\ of orders 1 to `order`, by whichever route the
class uses.

## Usage

``` r
cdf_tables(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- order:

  The highest order wanted, 1 to 4.

## Value

A list of length `order` of named derivative tables of \\F\\.

## Details

The conversion to the log scale needs every order below the one wanted,
not just the one wanted, because the moment-to-cumulant relation is a
sum over partitions of the multi-index and a partition into \\k\\ blocks
asks for \\k\\ lower-order ratios. Collecting them in one place keeps
the choice of route – exact sum for a discrete family, one stencil on
the analytic distribution function for a continuous one – in a single
statement.

## See also

[`cdf_scale_k`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md),
[`discrete_cdf_deriv_k`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md),
[`numerical_cdf_deriv_k`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md)
