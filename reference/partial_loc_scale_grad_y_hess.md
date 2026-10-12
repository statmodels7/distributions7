# Second-Order Mixed Derivatives With a Shape Parameter

The
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
method bodies for a family carrying shape parameters beyond its location
and scale. The three pairs in \\(\mu, \sigma)\\ come from
[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
in closed form, and every pair touching a shape parameter from one
central difference of the first-order component,
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the third-order body and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the fourth. In the fourth-order body a pair of shape parameters is
one stencil on the analytic
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
since the shape components of
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
are themselves a difference in the shape. Both delegate to
[`partial_theta2()`](https://statmodels7.github.io/distributions7/reference/partial_theta2.md),
which does the splicing.

## Usage

``` r
partial_loc_scale_grad_y_hess(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  ...
)

partial_loc_scale_hess_y_hess(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  ...
)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`, whose first two parameters
  are a location and a scale in that order and whose remaining
  parameters are shapes.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  One of `"parameter"` or `"link"`. The generic applies the link scale
  to the value that the method returns.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

Registered on the two pseudo-Huber parametrizations, the skew normal and
the skew t. The Student t's methods, documented on this page as well,
call compiled kernels of their own and do not use these bodies.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response, \\z =
(y-\mu)/\sigma\\ the standardized residual, and \\A\\ and \\B\\ the two
standardized quantities \\\sigma\ell^{(y)}\\ and
\\\sigma^2\ell^{(yy)}\\, whose derivatives are taken in \\z\\.

## See also

[`partial_theta2()`](https://statmodels7.github.io/distributions7/reference/partial_theta2.md),
which splices the two halves,
[`loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_y_hess.md)
for a family with no shape parameter, and
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the generics.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma = 1.2, alpha = 0.7)

g3 <- distrib_grad_y_hess(d, y, theta)
names(g3)
#> [1] "mu_mu"       "sigma_sigma" "alpha_alpha" "mu_sigma"    "mu_alpha"   
#> [6] "sigma_alpha"
vapply(g3, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma alpha_alpha    mu_sigma    mu_alpha sigma_alpha 
#>  0.03060209  4.59301048  1.09678158 -1.60525989  0.76774711 -1.77660202 

# Against a numerical Hessian of the response gradient. The three pairs in
# (mu, sigma) are closed and the three touching alpha are differenced.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2],
                                              alpha = v[3]))
numDeriv::hessian(f, c(0.3, 1.2, 0.7))
#>             [,1]      [,2]       [,3]
#> [1,]  0.03060209 -1.605260  0.7677471
#> [2,] -1.60525989  4.593010 -1.7766020
#> [3,]  0.76774711 -1.776602  1.0967816

# The closed three are exactly the location-scale block's.
cl <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
all(vapply(names(cl), function(k) identical(g3[[k]], cl[[k]]), logical(1)))
#> [1] TRUE
```
