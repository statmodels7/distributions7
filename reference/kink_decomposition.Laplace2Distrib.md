# The Rate-Parametrized Laplace's Kink

The log-density is \\\log(\lambda/2) - \lambda \lvert y - \mu \rvert\\,
so the composition is the absolute value at \\v = y - \mu\\ with
\\c(\theta) = -\lambda\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

## Value

A
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md)
with `phi = "abs"`.

## Details

The same kink as
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)'s,
carried across the map \\\lambda = 1/\sigma\\: a reparametrization moves
the coefficient in front and leaves the argument alone, so the order is
unchanged and the location is again the only non-smooth parameter.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md)
for the generic,
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
for the family.

## Examples

``` r
kink_decomposition(laplace2_distrib())
#> <kink_spec> c(theta) * phi(v),  phi = abs,  order 0
params_order(laplace2_distrib())
#>     mu lambda 
#>      0    Inf 
```
