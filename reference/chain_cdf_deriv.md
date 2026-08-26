# The Chain Rule on a Parent's CDF Derivatives

Carries the derivatives of a parent's distribution function onto a new
parametrization, given the map's partial derivatives. A family written
as a map of another therefore needs no cdf derivatives of its own: the
parent's closed forms and the map's partials, which the family already
declares for its likelihood derivatives, are between them enough.

## Usage

``` r
chain_cdf_deriv(parent, q, th_par, maps, new_params, order)
```

## Arguments

- parent:

  The distribution being reparametrized.

- q:

  A numeric vector of quantiles.

- th_par:

  The parent's parameters, evaluated at the new ones.

- maps:

  The map's keyed partial tables, in the form
  [`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
  produces: a missing key is an exact zero, so a map with many vanishing
  partials costs nothing for them.

- new_params:

  A character vector naming the new parameters, which names and orders
  the result.

- order:

  The derivative order, 1 or 2.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale: one per new parameter at order 1, and one per
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
component at order 2.

## The identity

The region of integration is independent of the parameters in both
parametrizations, so the ordinary chain rule applies to the composition
\\F(q; h(\psi))\\: \$\$\frac{\partial F}{\partial\psi_a} = \sum_k
F_k\\h^k_a, \qquad \frac{\partial^2 F}{\partial\psi_a \partial\psi_b} =
\sum\_{k,l} F\_{kl}\\h^k_a h^l_b + \sum_k F_k\\h^k\_{ab}.\$\$

## When it is worth taking

The result is exact whenever the parent's own cdf derivatives are, and
no better than they are otherwise: applied to a parent that differences
its own cdf, the chain adds an exact transformation to an approximate
quantity and buys nothing over differencing the new cdf directly.
Checking that is the caller's business, and
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md)
is where it is done.

## Notation

\\F\\ is the parent's distribution function, \\\psi\\ the new
parameters, \\h\\ the map from them to the parent's, \\F_k\\ and
\\F\_{kl}\\ the parent's cdf derivatives and \\h^k_a\\, \\h^k\_{ab}\\
the map's partials.

## See also

[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md),
the caller that gates it;
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
the same construction for likelihood derivatives;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).
