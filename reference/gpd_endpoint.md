# The Upper Endpoint of a Generalized Pareto

Returns \\-\sigma/\xi\\ where \\\xi \< 0\\ and `Inf` elsewhere: the
point beyond which a generalized Pareto puts no mass.

## Usage

``` r
gpd_endpoint(sigma, xi)
```

## Arguments

- sigma:

  A numeric vector of scales, strictly positive. Nothing is validated.

- xi:

  A numeric vector of shapes, of any sign.

## Value

A numeric vector of endpoints, of the length of the recycled inputs.
Entries are `Inf` wherever `xi >= 0`.

## Details

The endpoint depends on the parameters, which is the whole reason this
family needs care. For \\\xi \< 0\\ the support is \\\[0,
-\sigma/\xi\]\\ and both ends of the interval move as the parameters do,
so differentiating an expectation under the integral sign is not
automatically licensed. The consequence a user meets is that the
expected information exists only for \\\xi \> -1/2\\; see
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

[`ifelse()`](https://rdrr.io/r/base/ifelse.html) is correct here and is
deliberate: both arguments are of the length of the test after
recycling, so a vector of shapes gives a vector of endpoints.

## See also

[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
for the family and what the moving endpoint costs, and
[`distrib_expected_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GPDDistrib.md)
for where it stops existing.

## Examples

``` r
# Finite only for a negative shape.
distributions7:::gpd_endpoint(2, c(-0.4, -0.1, 0, 0.3))
#> [1]   5  20 Inf Inf

# The density really is zero at and beyond it.
d <- gpd_distrib()
distrib_pdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
#> [1] 0.04472136 0.00000000 0.00000000
```
