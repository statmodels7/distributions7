# CDF Derivatives of a Location-Scale Family

Closed-form derivatives of \\F\\ for a family that is location-scale in
its first two parameters.

## Usage

``` r
loc_scale_cdf_deriv(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, location first and scale second.

- order:

  The derivative order, 1 or 2.

## Value

A named list of derivative component vectors of \\F\\.

## Details

With \\z = (q - \mu)/\sigma\\ and \\\ell_y = \partial \log f/\partial
y\\, \$\$\partial F/\partial \mu = -f, \qquad \partial F/\partial \sigma
= -z f\$\$ \$\$\partial^2 F/\partial \mu^2 = f \ell_y, \qquad \partial^2
F/\partial \mu \partial \sigma = f (z \ell_y + 1/\sigma), \qquad
\partial^2 F/\partial \sigma^2 = f (z^2 \ell_y + 2z/\sigma).\$\$

A family therefore only has to declare that it is location-scale;
nothing is needed beyond its density and its response derivative, both
of which every distribution already provides. This covers the
censored-regression workhorses – Gaussian, logistic, Cauchy, Laplace.

## See also

[`loc_scale_grad_cdf`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md),
[`loc_scale_hess_cdf`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
