# Index Tuples Behind the Higher-Order Derivative Names

The multi-indices
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names, in exactly the same order: non-decreasing tuples of length
`order` over `seq_along(params)`, enumerated lexicographically.

## Usage

``` r
deriv_indices(params, order)
```

## Arguments

- params:

  A character vector of parameter names.

- order:

  A positive integer, the derivative order.

## Value

A list of integer vectors, each of length `order`, parallel to
`deriv_names(params, order)`.

## Details

This exists so that nothing has to recover an index tuple by splitting a
component name back apart. That is the obvious route and it is wrong: a
parameter whose own name contains an underscore makes
`"mu_log_scale_log_scale"` split into five pieces, and matching those
against `params` yields `NA`s. The failure is not subtle when it arrives
– the numerical third and fourth derivatives and the Bartlett expected
Hessian all stopped with an error – but it only arrives for a
user-defined distribution, since none of the fourteen shipped here has
an underscore in a parameter name. Generating the indices and the names
from the same enumeration cannot be fooled.

Note that this is **not** interchangeable with
[`deriv_index_list()`](https://statmodels7.github.io/distributions7/reference/deriv_index_list.md)
in `link_scale.R`: that one is ordered to match
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2, which puts the diagonal first, whereas this one is
lexicographic throughout to match
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).
At order 2 use
[`hess_pairs`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md).

## See also

[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
[`hess_pairs`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
