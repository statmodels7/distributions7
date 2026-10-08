# Third and Fourth CDF Derivatives of a Continuous Family

Returns the derivatives of \\F\\ of order three or four by
[`compiled_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_deriv_k.md)
for the families it covers, and by the stencil of
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md)
otherwise.

## Usage

``` r
continuous_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  A family named in the description.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The derivative order, 3 or 4.

## Value

A named list of numeric vectors, as
[`compiled_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_deriv_k.md).
