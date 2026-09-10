# The Laplace's Kink

The log-density is \\-\log(2\sigma) - \lvert y - \mu \rvert / \sigma\\,
so the composition is the absolute value at \\v = y - \mu\\ with
\\c(\theta) = -1/\sigma\\.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

## Value

A
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md)
with `phi = "abs"`.

## Details

Only the location moves the kink: \\\partial v/\partial \mu = -1\\ and
\\\partial v/\partial \sigma = 0\\, so
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
reports \\m\_\mu = 0\\ and \\m\_\sigma = \infty\\, which is what
`params_smooth` has recorded as `c(mu = FALSE, sigma = TRUE)` all along.

The scale is smooth even though it stands in front of the absolute
value, because a coefficient multiplying \\\phi\\ is differentiable
wherever \\\phi\\ is bounded; what breaks a derivative is the argument
crossing the origin, not the factor.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md)
for the generic,
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for the family.

## Examples

``` r
kink_decomposition(laplace_distrib())
#> <kink_spec> c(theta) * phi(v),  phi = abs,  order 0
params_order(laplace_distrib())
#>    mu sigma 
#>     0   Inf 
check_kink(laplace_distrib())
#> check_kink: laplace,  phi = abs,  order 0
#>   dv against a difference of v      [PASSED]  worst 4.55e-12
#>   the jump of the score             [PASSED]  worst 0.00e+00
#>   the smooth parameters             [PASSED]
#>       mu         measured            2   declared            2
#>       sigma      measured            0   declared            0
```
