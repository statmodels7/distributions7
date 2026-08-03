# Index Pairs Behind the Hessian Keys of a Multivariate Distribution

The positions in `distrib@params` of each unordered pair, in the order
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
uses.

## Usage

``` r
mv_hess_indices(distrib)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

## Value

A list of integer vectors of length 2.

## Details

[`hess_pairs`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
returns pairs of parameter NAMES, which is what a univariate method
wants when it looks a component up. A multivariate method needs the
positions instead, to tell a mean component from a matrix one, so the
names are matched back rather than the enumeration being written a
second time.
