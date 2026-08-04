# Affine (Location-Scale) Transformation

Transformer for \\Y = \text{loc} + \text{scale} \cdot X\\. Inverse \\X =
(Y - \text{loc})/\text{scale}\\, Jacobian \\\|J\| =
1/\|\text{scale}\|\\.

## Usage

``` r
affine_transform(loc = 0, scale = 1)
```

## Arguments

- loc:

  Numeric. The location shift. Defaults to 0.

- scale:

  Numeric. The scale multiplier (non-zero). Defaults to 1.

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gaussian_distrib(), affine_transform(loc = 1, scale = 2))
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.1994711
```
