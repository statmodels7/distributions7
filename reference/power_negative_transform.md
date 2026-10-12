# Even Power of a Negative Variable

Returns the transformer of \\Y = X^p\\ for an even integer \\p\\ and a
variable \\X \le 0\\, whose inverse is the negative root \\X =
-Y^{1/p}\\.
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
substitutes it for
[`power_transform()`](https://statmodels7.github.io/distributions7/reference/power_transform.md)
when the parent's support lies at or below zero. The map decreases for
\\p \> 0\\ and increases for \\p \< 0\\; the absolute Jacobian and its
logarithmic derivatives are those of the positive root.

## Usage

``` r
power_negative_transform(p)
```

## Arguments

- p:

  An even integer, the exponent.

## Value

A
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

[`power_transform()`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(truncated(gaussian1_distrib(), upper = 0),
                    power_transform(2))
d@bounds
#> [1]   0 Inf
distrib_pdf(d, 4, list(mu = -1, sigma = 1))
#> [1] 0.07189999
dnorm(-2, -1, 1) / pnorm(0, -1, 1) / 4
#> [1] 0.07189999
```
