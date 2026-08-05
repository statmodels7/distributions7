# The Partition Sum Itself

Carries the parent's derivatives into new coordinates, given the jets of
the map. Separated from
[`reparam_chain`](https://statmodels7.github.io/distributions7/reference/reparam_chain.md)
so that a family written in its own right, rather than obtained through
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
can use the same machinery instead of a second copy of it.

## Usage

``` r
chain_derivatives(parent, y, th_par, jt, new_params, order, expected = FALSE)
```

## Arguments

- parent:

  The distribution whose derivatives are being carried.

- y:

  The response.

- th_par:

  The parent's parameters, as plain numbers.

- jt:

  The jets of the map, as
  [`reparam_jets`](https://statmodels7.github.io/distributions7/reference/reparam_jets.md)
  returns them.

- new_params:

  The names of the new parameters.

- order:

  The derivative order, 1 to 4.

- expected:

  Logical; if `TRUE`, carries the expected derivatives.

## Value

A named list of component vectors.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
