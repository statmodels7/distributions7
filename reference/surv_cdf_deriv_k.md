# CDF Derivatives From an Exponential Survival Function

Returns \\\partial^I F\\ for every component of the requested order, for
a family whose survival function is the exponential of something
elementary. Given \\L = \log(1-F)\\ and a function evaluating its
partial derivatives, all four orders follow at once, so a family has
only to say what \\L\\ is.

## Usage

``` r
surv_cdf_deriv_k(distrib, q, theta, order, Lval, Lderiv, inside = NULL)
```

## Arguments

- distrib:

  An object inheriting from `distrib`. Its `params` name and order the
  components, and its lower bound supplies the default mask.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The derivative order, 1 to 4.

- Lval:

  The value of \\L\\ at `q`, a numeric vector.

- Lderiv:

  A function of a character vector of parameter names returning the
  corresponding partial derivative of \\L\\. An empty block is the
  zeroth order and is never asked for.

- inside:

  A logical vector saying which quantiles lie inside the support, or
  `NULL` (the default), which reads `q > distrib@bounds[1]`. A family
  whose support depends on a parameter, as the generalized Pareto's does
  at a negative shape, supplies its own; the fixed bounds cannot see it.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
and exactly zero wherever `inside` is `FALSE`.

## The identity

\\S = e^{L}\\ gives \\\partial^I S = S\\B_I\\, with \\B_I\\ the complete
Bell polynomial in the partials of \\L\\, and \\F = 1 - S\\ turns that
into \\\partial^I F = -S\\B_I\\. It is the same identity the
distribution wrappers use, read on the survival function instead of on
the density, and
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
runs the partition sum, so nothing here is transcribed from an
expansion.

## Two things the arithmetic needs

The survival function is evaluated as `exp(Lval)` and not as `1 - F`,
which keeps the far tail from canceling. And below the support \\F\\ is
identically zero and so is every derivative, while \\L\\ is still finite
there and would otherwise give a survival above one; the `inside` mask
is what suppresses that.

## Notation

\\F\\ is the distribution function, \\S = 1 - F\\ the survival function,
\\L = \log S\\ and \\B_I\\ the complete Bell polynomial.

## See also

[`register_surv_cdf()`](https://statmodels7.github.io/distributions7/reference/register_surv_cdf.md),
which turns a pieces function into the four methods;
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
for the partition sum;
[`gpd_surv_pieces()`](https://statmodels7.github.io/distributions7/reference/gpd_surv_pieces.md)
for the most involved of the three families.
