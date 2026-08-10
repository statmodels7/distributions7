# The Location and Scale Block of a Second-Order Mixed Derivative

The three components in \\(\mu, \sigma)\\ of
\\\partial^2\ell^{(y)}/\partial\theta^2\\ or of
\\\partial^2\ell^{(yy)}/\partial\theta^2\\, for a family whose response
enters only through \\z = (y-\mu)/\sigma\\.

## Usage

``` r
loc_scale_theta2_block(distrib, y, theta, order = 1L)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, location first and scale second.

- order:

  `1` for the block of
  [`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` for that of
  [`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).

## Value

A named list of three components, keyed `mu_mu`, `sigma_sigma`,
`mu_sigma` under the family's own parameter names.

## Details

The quantities are the family's own response derivatives scaled by a
power of \\\sigma\\: writing \\A = \sigma\ell^{(y)}\\, its derivatives
in \\z\\ are \\\sigma^2\ell^{(yy)}\\ and \\\sigma^3\ell^{(yyy)}\\, and
for \\B = \sigma^2\ell^{(yy)}\\ they are \\\sigma^3\ell^{(yyy)}\\ and
\\\sigma^4\ell^{(yyyy)}\\. Nothing new is differentiated.

## See also

[`distrib_grad_y_hess`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
