# Mixed Derivatives by One Stencil on the Highest Analytical Quantity

Computes the derivatives of the log-density of order `ay` in the
response and `ctheta` in the parameters, for a continuous family that
has no method of its own for them. Among the quantities the family
implements itself (the log-density, its response derivatives, its score
and Hessian in the parameters, and the mixed derivatives of lower
order), the one that leaves the smallest order to difference is taken,
and each component is a single tensor-product central difference of it,
never a difference of a difference.

## Usage

``` r
mixed_tensor_derivatives(distrib, y, theta, ay, ctheta)
```

## Arguments

- distrib:

  A `continuous_distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- ay:

  The order in the response, 1 to 3.

- ctheta:

  The order in the parameters, 1 or 2.

## Value

A named list: one vector per parameter, keyed by `distrib@params`, where
`ctheta` is 1; one vector per unordered pair, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
where it is 2.

## Details

A candidate quantity of order \\(a, c)\\ in the response and the
parameters can serve where \\a \le\\ `ay` and \\c \le\\ `ctheta`, and it
leaves \\r = (\mathrm{ay} - a) + (\mathrm{ctheta} - c)\\ orders to the
stencil. The log-density always serves; another quantity serves only
where the family's class registers its method, as
[`owns_method()`](https://statmodels7.github.io/distributions7/reference/owns_method.md)
reads it, so that a numerical fallback is never differenced again. The
candidate with the smallest \\r\\ is used, and among equal ones the one
with more derivatives in the parameters.

The stencil is the tensor product of one-dimensional central stencils of
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
at second-order accuracy, one in the response and one in each parameter
that still has to be differenced, with the parameters that occur most
often taken by the analytical quantity first, as
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md)
does. The step is \\h = \varepsilon^{1/(r+2)}\max(1, \|x\|)\\ in every
direction, the response step clamped inside the support by
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
and a parameter's step inside its bounds. Values shared by several
components are computed once.

For
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which carries the parent's score and Hessian in the parameters and no
response derivative, the fourth-order \\\partial^4\ell/\partial
y^2\\\partial\theta^2\\ is a second difference in the response of the
analytical Hessian: \\r = 2\\, where the previous fallback differenced a
difference in the same parameter.

## See also

[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md)
for the derivatives in the parameters alone,
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
whose continuous fallbacks call this.

## Examples

``` r
# transformation() has no response derivative of its own: the mixed
# fourth order of a log-gamma is a stencil on its analytical Hessian.
d <- fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)
b <- c(-0.4, 0.3)
k <- 1 / 0.26
got <- distributions7:::mixed_tensor_derivatives(d, b, list(sigma2 = 0.26), 2, 2)
#> Error in vapply(usable, function(b) (ay - b$a) + (ctheta - b$c), integer(1)): values must be type 'integer',
#>  but FUN(X[[1]]) result is type 'double'
rbind(got$sigma2_sigma2, -2 * k^3 * exp(b))
#> Error: object 'got' not found
```
