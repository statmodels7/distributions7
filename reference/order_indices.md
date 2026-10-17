# Multi-Indices of a Given Order, as Parameter Names

Returns the multi-indices of a given order as parameter names, in the
order in which
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
lists them.

## Usage

``` r
order_indices(params, order)
```

## Arguments

- params:

  A character vector of parameter names.

- order:

  The derivative order.

## Value

A list of character vectors, each of length `order`.

## Details

This function is a thin wrapper on
[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md).
It is deliberately not
[`deriv_index_list()`](https://statmodels7.github.io/distributions7/reference/deriv_index_list.md)
from `link_scale.R`, whose order-2 case is ordered for
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
– diagonal first – while
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
is lexicographic; pairing those would silently attach the name
`"mu_sigma"` to the index `(sigma, sigma)`. The orders used here are 3,
4 and 5, where the two orderings agree; the function generates the
indices locally, which avoids the mismatch and does not depend on order
2 never being requested.

## See also

[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md),
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
