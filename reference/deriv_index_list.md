# Index Tuples Matching the Package's Component Naming

Canonical (non-decreasing) index tuples for a derivative order, in
exactly the output order of the corresponding name helper: parameter
order at order 1,
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2,
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
from order 3 up.

## Usage

``` r
deriv_index_list(p, order)
```

## Arguments

- p:

  The number of parameters.

- order:

  The derivative order, 1 to 4.

## Value

A list of integer vectors, each of length `order`.

## Details

The order-2 case is the one to be careful with.
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
lists the diagonal first and the off-diagonal afterwards, whereas
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
is lexicographic throughout; pairing this helper with
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
would therefore label `"mu_sigma"` with the tuple `(sigma, sigma)`.
Orders 3 and 4 agree between the two conventions, so the mismatch is
invisible until someone reuses the helper at order 2. Use
[`deriv_indices`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
when the names come from
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
and
[`hess_pairs`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
when they come from
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`deriv_indices`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md),
[`hess_pairs`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
