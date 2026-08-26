# CDF Derivatives of a Location-Scale Family

Closed-form derivatives of \\F\\ for a family that is location-scale in
its first two parameters. Writing \\F(q) = G(z)\\ with \\z =
(q-\mu)/\sigma\\ and differentiating the composition puts everything in
terms of the density and its response derivative, both of which every
distribution already supplies, so a family has only to declare that it
is location-scale.

## Usage

``` r
loc_scale_cdf_deriv(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `distrib` whose first two parameters are a
  location and a scale. Nothing checks that; a family that is not
  location-scale and registers this body gets wrong numbers.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, the location first and the scale second.
  Any further parameters are read by
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
  and
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  but get no component here.

- order:

  The derivative order, 1 or 2.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale: two components at order 1, and the three of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2.

## Details

With \\\ell_y = \partial \log f/\partial y\\, \$\$\frac{\partial
F}{\partial \mu} = -f, \qquad \frac{\partial F}{\partial \sigma} = -z
f,\$\$ \$\$\frac{\partial^2 F}{\partial \mu^2} = f\\\ell_y, \qquad
\frac{\partial^2 F}{\partial \mu \partial \sigma} = f\left(z\\\ell_y +
\frac{1}{\sigma}\right), \qquad \frac{\partial^2 F}{\partial \sigma^2} =
f\left(z^2 \ell_y + \frac{2z}{\sigma}\right).\$\$

This covers the censored-regression workhorses: the Gaussian, the
logistic, the Cauchy and the Laplace all register these bodies
unchanged, and the Student t and the pseudo-Huber use them for their
first two components and difference the shape.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\\ell_y = \partial\log
f/\partial y\\ its response derivative.

## See also

[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
and
[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md),
the two bodies families register;
[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md),
the variant for a family with a shape parameter as well;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).
