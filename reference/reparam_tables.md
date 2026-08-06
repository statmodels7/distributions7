# The Map Derivatives as Keyed Tables

Returns, for every parent parameter, the partial derivatives of the map
component with respect to the new parameters, keyed by the sorted tuple
of new-parameter positions ("1", "1,2", "2,2,3,3", ...). A missing key
is an exact zero.

## Usage

``` r
reparam_tables(distrib, theta)
```

## Arguments

- distrib:

  A reparametrized distribution.

- theta:

  A named list of the new parameters, already aligned.

## Value

A list over parent parameters of keyed partial tables.

## Details

When the family supplies `map_derivs`, the tables are its hand-written
closed forms; the shipped second parametrizations all do, and the
formulas live in
[`reparam_map_derivs`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md).
Otherwise each needed partial comes from one finite-difference stencil
of
[`fd_derivative`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
applied to the analytic map – a single stencil per order, never a chain
of differences, at the accuracy that construction carries (about 1e-8 at
first order, fading with the order). Exact tables are therefore the
recommendation for any family fitted in earnest.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
