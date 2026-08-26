# One Stencil Per Map Partial

The numerical fallback behind
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md):
every partial of the map is one central stencil of
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
in each direction, applied to the analytic map, iterated across
**distinct** directions. A cross-variable composition of two first
differences is one mixed stencil and not the same-variable nesting the
toolkit forbids.

## Usage

``` r
reparam_stencil_derivs(map, params, parent_params)
```

## Arguments

- map:

  The map function, of a named list of the new parameters returning a
  named list of the parent's.

- params:

  A character vector naming the new parameters.

- parent_params:

  A character vector naming the parent's, which fixes the order of the
  returned list.

## Value

A function of the new parameters, usable as an object's
`reparam_derivs`, returning the keyed tables
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
describes.

## Details

The accuracy is about \\10^{-8}\\ at first order and fades with the
order, which is why a family fitted in earnest supplies `map_derivs`
instead. What the fallback buys is that a user's own reparametrization
works immediately, with no derivatives written: the same bargain
[`distrib_grad_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.continuous_distrib.md)
offers a density-only family.

## See also

[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md),
which calls the result;
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
for the stencil;
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
for the exact alternative.
