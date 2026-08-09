# CDF Derivatives of a Discrete Distribution, at Any Order

The general form of
[`discrete_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md):
evaluates \\d^I F(q) = \sum\_{y \le q} f(y)\\(d^I f/f)(y)\\ for any
order up to four.

## Usage

``` r
discrete_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

## Value

A named list of derivative components of \\F\\.

## Details

The quantity summed is the complete Bell polynomial in the
log-derivatives, which
[`bell_f_ratio`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
computes, so the order enters only through how many of the family's
derivative tables are fetched. At orders one and two this reproduces the
written-out \\f g\\ and \\f(h + g_i g_j)\\ of
[`discrete_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md).

As there, a test must not check this against the partial-expectation
sum, which is the same sum twice; finite differences of the cdf are the
independent reference.

## See also

[`discrete_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md),
[`bell_f_ratio`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
