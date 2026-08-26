# The Location and Scale Block of a Second-Order Mixed Derivative

Computes the three components in \\(\mu, \sigma)\\ of
\\\partial^2\ell^{(y)}/\partial\theta^2\\ or of
\\\partial^2\ell^{(yy)}/\partial\theta^2\\, for a family whose response
enters only through \\z = (y-\mu)/\sigma\\. It is the arithmetic behind
every location-scale method of
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
and it derives nothing new.

## Usage

``` r
loc_scale_theta2_block(distrib, y, theta, order = 1L)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`, whose first two parameters
  are a location and a scale in that order.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, location first and scale second.

- order:

  `1` for the block of
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` for that of
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).
  Default `1`.

## Value

A named list of three numeric vectors, keyed `mu_mu`, `sigma_sigma` and
`mu_sigma` under the family's own parameter names.

## The identity

With \\A = \sigma\ell^{(y)}\\ and \\B = \sigma^2\ell^{(yy)}\\ functions
of \\z\\ alone, and \\\partial z/\partial\mu = -1/\sigma\\, \\\partial
z/\partial\sigma = -z/\sigma\\, the chain rule gives
\$\$\frac{\partial^2\ell^{(y)}}{\partial\mu^2} = \frac{A''}{\sigma^3},
\qquad \frac{\partial^2\ell^{(y)}}{\partial\mu\\\partial\sigma} =
\frac{zA'' + 2A'}{\sigma^3}, \qquad
\frac{\partial^2\ell^{(y)}}{\partial\sigma^2} = \frac{z^2A'' + 4zA' +
2A}{\sigma^3},\$\$ and at the next order in the response
\$\$\frac{\partial^2\ell^{(yy)}}{\partial\mu^2} = \frac{B''}{\sigma^4},
\qquad \frac{\partial^2\ell^{(yy)}}{\partial\mu\\\partial\sigma} =
\frac{zB'' + 3B'}{\sigma^4}, \qquad
\frac{\partial^2\ell^{(yy)}}{\partial\sigma^2} = \frac{z^2B'' + 6zB' +
6B}{\sigma^4}.\$\$

## Why nothing new is differentiated

\\A\\, \\A'\\, \\A''\\ and \\B\\, \\B'\\, \\B''\\ are the family's OWN
response derivatives times a power of \\\sigma\\: \\A =
\sigma\ell^{(y)}\\, \\A' = B = \sigma^2\ell^{(yy)}\\, \\A'' = B' =
\sigma^3\ell^{(yyy)}\\ and \\B'' = \sigma^4\ell^{(yyyy)}\\. A family
that already carries
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
and
[`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
therefore gets both orders in the parameters for free, which is the
bargain
[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
already takes at first order.

A shape parameter beyond the two is NOT covered and falls back; see
[`partial_theta2()`](https://statmodels7.github.io/distributions7/reference/partial_theta2.md),
which splices these three components into the differenced ones.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response, \\z =
(y-\mu)/\sigma\\ the standardized residual, and \\A\\ and \\B\\ the two
standardized quantities \\\sigma\ell^{(y)}\\ and
\\\sigma^2\ell^{(yy)}\\, whose derivatives are taken in \\z\\.

## See also

[`loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_y_hess.md)
and
[`partial_loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md),
the two method bodies built on it,
[`loc_scale_cross2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_block.md)
for the first order in \\\theta\\, and
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma = 1.2)

blk <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
blk$mu_mu
#> [1] -0.09632533  0.00000000  0.10125089

# The identity written out: A'' / sigma^3, with A'' the family's own third
# response derivative scaled.
s <- 1.2
A2 <- s^3 * distrib_deriv3_y(d, y, theta)
A2 / s^3
#> [1] -0.09632533  0.00000000  0.10125089

# And the scale pair, which carries all three of A, A' and A''.
z <- (y - 0.3) / s
A <- s * distrib_grad_y(d, y, theta)
A1 <- s^2 * distrib_hess_y(d, y, theta)
c(reported = blk$sigma_sigma[1],
  formula = ((z^2 * A2 + 4 * z * A1 + 2 * A) / s^3)[1])
#> reported  formula 
#> 1.203953 1.203953 

# Order 2 reads one derivative further and is the same shape.
distributions7:::loc_scale_theta2_block(d, y, theta, 2L)$mu_mu
#> [1] 0.05438174 0.12056327 0.04414217
(s^4 * distrib_deriv4_y(d, y, theta)) / s^4
#> [1] 0.05438174 0.12056327 0.04414217
```
