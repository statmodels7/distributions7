# Derivatives by One Stencil on the Highest Analytical Order

Computes the derivative components of order `order` of the log-density
from the highest order the family implements itself: the log-density
(`base = 0`), the score (`base = 1`), the Hessian (`base = 2`), the
third order (`base = 3`) or the fourth (`base = 4`). Each component is a
single tensor-product central difference of order `order - base`, never
a difference of a difference.

## Usage

``` r
tensor_derivatives(
  distrib,
  y,
  theta,
  order,
  base,
  h_rel = NULL,
  skip = NULL,
  scale = c("parameter", "link"),
  accuracy = 4L
)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- order:

  The order of the derivatives, a whole number above `base`.

- base:

  The highest order the family implements itself, 0 to 4.

- h_rel:

  `NULL`, the default, for the step above, or a relative step that
  replaces \\\varepsilon^{1/(r+a)}\\.

- skip:

  Character vector of component names left `NULL`, or `NULL`.

- scale:

  `"parameter"` or `"link"`, the scale of the coordinates and of the
  result.

- accuracy:

  The order of accuracy of the one-dimensional stencils, passed to
  [`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html).

## Value

A named list of component vectors, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

The computation is
[`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html)
over the coordinates \\\theta_1, \dots, \theta_p\\, with the log-density
as the function and the derivatives of order `base` declared as the
known quantity. A component is a multiset of parameters: \\(\mu, \mu,
\phi)\\ is \\\partial^3\ell / \partial\mu^2\partial\phi\\. `base` of its
indices select the analytical component that is differenced, taken from
the parameters that occur most often, and the remaining multiplicities
are the orders of the one-dimensional central stencils whose tensor
product is applied to it. Each parameter is differentiated only once,
however many times it occurs.

The step in \\\theta_j\\ is \\h_j = \varepsilon^{1/(r+a)}\max(1,
\|\theta_j\|)\\, with \\r = \\ `order - base` the order of the
difference and \\a\\ the accuracy, reduced so that every node stays
inside the parameter's bounds. On the link scale the coordinates are the
linear predictors \\\eta_j\\, the bounds are those of the links, and the
known quantity is read on the link scale, so the result is the
derivative of a link-scale quantity.

## See also

[`analytic_order()`](https://statmodels7.github.io/distributions7/reference/analytic_order.md)
for how `base` is found;
[`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md),
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
and
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
which call this.
