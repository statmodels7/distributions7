# Inverse-Link Derivatives for Every Parameter

Returns the derivatives \\h', h'', \dots, h^{(6)}\\ of each parameter's
inverse link, evaluated at \\\eta = g(\theta)\\, up to `order`.

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

  The highest derivative order needed, 1 to 6.

## Value

A list with one element per parameter; each element is a list whose
\\k\\-th element is \\h^{(k)}\\ evaluated at that parameter's \\\eta\\.

## Details

Each link is evaluated by one call to
[`linkfunctions7::linkinv_derivs()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv_derivs.html),
which returns its derivatives of orders one to `order` together. See
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
for what the derivatives are then used for.

## See also

[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
