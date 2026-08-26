# The Map Derivatives as Keyed Tables

Returns, for every parent parameter, the partial derivatives of that
component of the map with respect to the new parameters, keyed by the
sorted tuple of new-parameter positions: `"1"`, `"1,2"`, `"2,2,3,3"` and
so on. A missing key is an exact zero, so a map with many vanishing
partials costs nothing for them.

## Usage

``` r
reparam_tables(distrib, theta)
```

## Arguments

- distrib:

  A reparametrized distribution.

- theta:

  A named list of the new parameters, already aligned by
  [`reparam_theta()`](https://statmodels7.github.io/distributions7/reference/reparam_theta.md)'s
  caller. Only the first `n_params` components are read.

## Value

A list over the parent's parameters, each element a keyed list of that
component's partial derivatives.

## Where the tables come from

When the family supplies `map_derivs`, they are its hand-written closed
forms; the shipped second parametrizations all do, and the formulas live
in
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md).
Otherwise each needed partial comes from one finite-difference stencil
of
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
applied to the analytic map, a single stencil per order and never a
chain of differences, at the accuracy that construction carries: about
\\10^{-8}\\ at first order, fading with the order. Exact tables are
therefore the recommendation for any family fitted in earnest.

## See also

[`reparam_stencil_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_stencil_derivs.md)
for the numerical route;
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
for the hand-written tables;
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
the consumer.
