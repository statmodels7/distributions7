# Mixed Derivatives by One Stencil on the Highest Analytical Quantity

Computes the derivatives of the log-density of order `ay` in the
response and `ctheta` in the parameters, for a continuous family that
has no method of its own for them. Among the quantities that the family
implements itself (the log-density, its response derivatives, its score
and Hessian in the parameters, and the mixed derivatives of lower
order), the one that leaves the smallest order to difference is taken,
and each component is a single tensor-product central difference of it,
never a difference of a difference.

## Usage

``` r
mixed_tensor_derivatives(distrib, y, theta, ay, ctheta, accuracy = 4L)
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

- accuracy:

  The order of accuracy of the one-dimensional stencils, passed to
  [`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html).

## Value

A named list: one vector per parameter, keyed by `distrib@params`, where
`ctheta` is 1; one vector per unordered pair, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
where it is 2.

## Details

The computation is
[`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html)
over the coordinates \\(y, \theta_1, \dots, \theta_p)\\. A candidate
quantity of order \\(a, c)\\ in the response and the parameters is
declared to it where \\a \le\\ `ay`, \\c \le\\ `ctheta`, \\(a, c)\\ is
not the target itself, and the family's class registers its method, as
[`owns_method()`](https://statmodels7.github.io/distributions7/reference/owns_method.md)
reads it, so that a numerical fallback is never differenced again. The
log-density is always available. The candidate with the smallest
residual order \\r = (\mathrm{ay} - a) + (\mathrm{ctheta} - c)\\ is
used, and among equal ones the candidate with more derivatives in the
parameters, which declares more multi-indices, is used; the candidates
are declared in decreasing \\c\\, so that the order of declaration gives
the same preference when there is a single parameter.

The step is \\h = \varepsilon^{1/(r+a)}\max(1, \|x\|)\\ in every
direction, with \\a\\ the accuracy, reduced so that every node stays
inside the support in the response and inside the bounds of each
parameter. Values shared by several components are computed once.

For
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which carries the parent's score and Hessian in the parameters and no
response derivative, the fourth-order \\\partial^4\ell/\partial
y^2\\\partial\theta^2\\ is a second difference in the response of the
analytical Hessian, so \\r = 2\\.

## See also

[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md)
for the derivatives in the parameters only,
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
got <- distributions7:::mixed_tensor_derivatives(d, b, list(sigma2 = 0.26), 2L, 2L)
rbind(got$sigma2_sigma2, -2 * k^3 * exp(b))
#>           [,1]      [,2]
#> [1,] -76.27675 -153.6025
#> [2,] -76.27675 -153.6025
```
