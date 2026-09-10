# The Elastic Net's Kink

The log-density is \\-a\lvert y-\mu\rvert - c(y-\mu)^2/2 - \log Z\\ with
\\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\, so the composition
is the absolute value at \\v = y - \mu\\ with \\c(\theta) =
-\lambda\alpha\\.

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

The SIZE of the kink is \\2\lambda\alpha\\ and goes to zero with
\\\alpha\\, which is why a detector reading the second Bartlett identity
finds this family in one sweep and not in another: at a small \\\alpha\\
the family is nearly Gaussian and the missing curvature disappears into
the scale of the rest of the matrix. The kink is there at every \\\alpha
\> 0\\ all the same, and the declaration says so where a measurement of
the curvature cannot.

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
