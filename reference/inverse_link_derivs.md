# Inverse-Link Derivatives for Every Parameter

The derivatives \\h', h'', h''', h''''\\ of each parameter's inverse
link, evaluated at \\\eta = g(\theta)\\, up to `order`.

## Usage

``` r
inverse_link_derivs(distrib, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameters, on the natural scale.

- order:

  The highest derivative order needed, 1 to 4.

## Value

A list with one element per parameter, each a list whose \\k\\-th
element is \\h^{(k)}\\ evaluated at that parameter's \\\eta\\.

## Details

This is the hot path of the link scale, so it is written against the
order-specific generics of linkfunctions7 rather than the convenience
router; see the comment in the body for why, and
[`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
for what the derivatives are then used for.

## See also

[`to_link_scale`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
