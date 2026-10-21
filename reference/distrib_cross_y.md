# Mixed Response-Parameter Derivatives of the Log-Density

Computes the mixed second derivatives \\\partial^2 \ell / \partial y\\
\partial \theta_i\\, with one component per parameter, and each
component is a vector along `y`. Together with
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
this completes the joint Hessian of the log-density in \\(y, \theta)\\:
those two are the diagonal blocks and this is the off-diagonal one.

The function is defined only for continuous families. A discrete family
has no derivative in \\y\\ to cross, so the package registers no method
for it and the call signals an error.

## Usage

``` r
distrib_cross_y(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  An object inheriting from `continuous_distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list, or a named numeric vector, of parameters. Each component
  must have length 1 or `length(y)`; a component of length 1 is
  recycled. The generic reorders `theta` by name and validates it
  against `params_bounds`, which are treated as open, before dispatch.

- scale:

  Either `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic applies it after dispatch, so a method always returns the
  parameter scale and never reads this argument.

- ...:

  Passed to the method. The shipped methods do not read it.

## Value

A named list with one numeric vector per parameter, keyed and ordered by
`distrib@params`. Each vector has length `length(y)`.

## What consumes it

A penalty in penalties7 is a negative log-density read at the
coefficients, so estimating coefficients and hyperparameters together
needs exactly this block, and so does the gradient of a profiled
objective through the implicit function theorem.

## The link scale

The component for \\\eta_i\\ is the parameter-scale component multiplied
by \\h_i'(\eta_i)\\, and nothing else. The response derivative is
untouched by a reparametrization of \\\theta\\, so only the first-order
diagonal chain rule enters, exactly as for
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md):
on a Gaussian whose scale carries a log link the `sigma` component is
multiplied by \\\sigma\\ and the `mu` component by 1.

## Where the closed forms are

Every continuous family in the package registers one. The fallback in
[`distrib_cross_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.continuous_distrib.md)
therefore serves families defined outside the package. It applies one
tensor stencil to the highest quantity that the family implements
itself, through
[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md).

## See also

[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
the two diagonal blocks;
[`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)
for the fallback;
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the third-order block.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
th <- list(mu = 0.3, sigma = 1.4)
distrib_cross_y(d, y, th)
#> $mu
#> [1] 0.5102041 0.5102041 0.5102041
#> 
#> $sigma
#> [1] -0.9475219 -0.2186589  1.2390671
#> 

# The mu component is 1/sigma^2 at every observation, and the sigma
# component is 2r/sigma^3, so it changes sign with the residual.
all.equal(distrib_cross_y(d, y, th)$sigma, 2 * (y - 0.3) / 1.4^3)
#> [1] TRUE

# On the link scale only the diagonal chain-rule factor enters: sigma
# carries a log link, so its component is multiplied by sigma itself.
distrib_cross_y(d, y, th, scale = "link")$sigma /
  distrib_cross_y(d, y, th)$sigma
#> [1] 1.4 1.4 1.4
```
