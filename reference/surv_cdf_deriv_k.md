# CDF Derivatives From an Exponential Survival Function

Returns \\\partial^{I}F\\ for every component of the requested order,
given \\L = \log(1-F)\\ and a function evaluating its partial
derivatives.

## Usage

``` r
surv_cdf_deriv_k(distrib, q, theta, order, Lval, Lderiv, inside = NULL)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

- Lval:

  The value of \\L\\ at `q`.

- Lderiv:

  A function of a character vector of parameter names returning the
  corresponding partial derivative of \\L\\.

## Value

A named list of derivative components of \\F\\.

## Details

\\S = e^{L}\\ gives \\\partial^{I}S = S\\B\_{I}\\, with \\B\_{I}\\ the
complete Bell polynomial in the partials of \\L\\, and \\F = 1 - S\\
turns that into \\\partial^{I}F = -S\\B\_{I}\\. The survival function is
evaluated as `exp(L)` rather than as `1 - F`, which keeps the far tail
from cancelling.

## See also

[`bell_f_ratio`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
