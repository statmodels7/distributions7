# The Elastic Net's Kink

The log-density is \\-a\lvert y-\mu\rvert - c(y-\mu)^2/2 - \log Z\\ with
\\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\, so the composition
is the absolute value at \\v = y - \mu\\ with coefficient
\\-\lambda\alpha = -a\\.

## Arguments

- distrib:

  An `EnetDistrib` object.

## Value

A
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md)
with `phi = "abs"`.

## Details

The quadratic half and the normalizing constant are smooth, so the
location is again the only non-smooth parameter, and \\m\_\mu = 0\\ for
every \\\alpha \> 0\\.

The size of the kink is \\2\lambda\alpha\\ and goes to zero with
\\\alpha\\. A detector based on the second Bartlett identity therefore
sees the missing curvature less clearly at small \\\alpha\\, where the
family is nearly Gaussian. The kink exists at every \\\alpha \> 0\\, and
the declaration records it independently of any measurement of the
curvature.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md)
for the generic,
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
for the family.

## Examples

``` r
kink_decomposition(enet_distrib())
#> <kink_spec> c(theta) * phi(v),  phi = abs,  order 0
params_order(enet_distrib())
#>     mu lambda  alpha 
#>      0    Inf    Inf 
```
