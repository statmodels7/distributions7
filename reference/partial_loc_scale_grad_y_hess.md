# Second-Order Mixed Derivatives With a Shape Parameter

The
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
method bodies for a family carrying shape parameters beyond its location
and scale. The three pairs in \\(\mu, \sigma)\\ come from
[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
in closed form, and every pair touching a shape parameter from one
central difference of the analytic first-order component. Both delegate
to
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

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

Registered on the Student t, the pseudo-Huber, the skew normal and the
skew t. Measured against a numerical Hessian of the response derivative,
the closed pairs agree to about \\10^{-11}\\ and the differenced ones to
between \\10^{-7}\\ and \\10^{-6}\\, which is one stencil's accuracy.

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
d <- student_t1_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma = 1.2, nu = 6)

g3 <- distrib_grad_y_hess(d, y, theta)
names(g3)
#> [1] "mu_mu"       "sigma_sigma" "nu_nu"       "mu_sigma"    "mu_nu"      
#> [6] "sigma_nu"   
vapply(g3, function(z) z[1], numeric(1))
#>        mu_mu  sigma_sigma        nu_nu     mu_sigma        mu_nu     sigma_nu 
#> -0.389443959  2.336663756  0.001414536 -0.634612372  0.018751535 -0.011680640 

# Against a numerical Hessian of the response gradient. The three pairs in
# (mu, sigma) are closed and the three touching nu are differenced.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2],
                                              nu = v[3]))
numDeriv::hessian(f, c(0.3, 1.2, 6))
#>             [,1]        [,2]         [,3]
#> [1,] -0.38944396 -0.63461237  0.018751535
#> [2,] -0.63461237  2.33666376 -0.011680640
#> [3,]  0.01875153 -0.01168064  0.001414536

# The closed three are exactly the location-scale block's.
cl <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
all(vapply(names(cl), function(k) identical(g3[[k]], cl[[k]]), logical(1)))
#> [1] TRUE
```
