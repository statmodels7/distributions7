# One Stencil Per Map Partial

The numerical fallback behind
[`reparam_tables`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md):
every partial of the map is one central stencil of
[`fd_derivative`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
in each direction, applied to the analytic map, iterated across DISTINCT
directions (a cross-variable composition, which is not the forbidden
same-variable nesting).

## Usage

``` r
reparam_stencil_derivs(map, params, parent_params)
```

## Arguments

- map:

  The map function.

- params:

  The new parameter names.

- parent_params:

  The parent parameter names.

## Value

A function usable as `reparam_derivs`.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
