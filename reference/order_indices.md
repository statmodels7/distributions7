# Multi-Indices of a Given Order, as Parameter Names

The multi-indices of a given order, in the order
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
lists them, expressed as parameter names.

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

A thin wrapper on
[`deriv_indices`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md).
It is deliberately not
[`deriv_index_list()`](https://statmodels7.github.io/distributions7/reference/deriv_index_list.md)
from `link_scale.R`, whose order-2 case is ordered for
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
– diagonal first – while
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
is lexicographic; pairing those would silently attach the name
`"mu_sigma"` to the index `(sigma, sigma)`. The orders actually
registered here are 3 and 4, where the two agree, but a mismatch that
only bites when someone reuses the helper is the kind worth removing
rather than commenting on.

## See also

[`deriv_indices`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md),
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
