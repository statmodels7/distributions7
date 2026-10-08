# Derivatives by One Stencil on the Highest Analytical Order

Computes the derivative components of order `order` of the log-density
from the highest order the family implements itself: the log-density
(`base = 0`), the score (`base = 1`), the Hessian (`base = 2`) or the
third order (`base = 3`). Each component is a single central difference
of order `order - base`, never a difference of a difference.

## Usage

``` r
tensor_derivatives(distrib, y, theta, order, base, h_rel = NULL, skip = NULL)
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

  The highest order the family implements itself, 0 to 3.

- h_rel:

  `NULL`, the default, for the step above, or a relative step that
  replaces \\\varepsilon^{1/(r+2)}\\.

- skip:

  Character vector of component names left `NULL`, or `NULL`.

## Value

A named list of component vectors, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

A component is a multiset of parameters: \\(\mu, \mu, \phi)\\ is
\\\partial^3\ell / \partial\mu^2\partial\phi\\. `base` of its indices
select the analytical component that is differenced, taken from the
parameters that occur most often so that no direction carries a
high-order stencil while another carries none; the remaining
multiplicities \\r_1, \dots, r_p\\ are the orders of the difference in
each direction. The stencil is the tensor product of one-dimensional
central stencils, the one for direction \\j\\ being
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
of order \\r_j\\ at second-order accuracy, so the component is
\$\$\sum\_{k_1, \dots, k_p} \Bigl(\prod_j w^{(r_j)}\_{k_j} /
h_j^{r_j}\Bigr) f(\theta_1 + k_1 h_1, \dots, \theta_p + k_p h_p),\$\$ a
single linear combination of values of the base component \\f\\. It
differentiates each parameter only once, however many times that
parameter occurs, which is what separates it from differencing a Hessian
that is itself a difference.

The step is \\h_j = \varepsilon^{1/(r+2)}\max(1, \|\theta_j\|)\\, with
\\r = \sum_j r_j\\ the total order of the difference: rounding grows as
\\\varepsilon/h^r\\ and the truncation of a second-order stencil as
\\h^2\\, and the two balance there. The step is clamped so that the
widest stencil of that order stays inside the parameter's bounds. Every
component of one order uses the same steps, so a point shared by several
components is evaluated once.

## See also

[`analytic_order()`](https://statmodels7.github.io/distributions7/reference/analytic_order.md)
for how `base` is found;
[`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md)
and
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
which call this for a family without its own Hessian.
