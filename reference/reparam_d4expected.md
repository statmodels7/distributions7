# The Fourth Derivative of the Expected Information of a Reparametrized Family

Reads the parent's expected information and its first four derivatives
at the mapped parameters and the map's tables to the fifth order, and
passes them to
[`dexpected_chain4()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain4.md).

## Usage

``` r
reparam_d4expected(distrib, y, theta, threads = 1L)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response, read for its length.

- theta:

  A named list of the new parameters.

- threads:

  The thread count passed to the parent's methods.

## Value

A named list on the parameter scale, keyed as
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

## Details

The tables of the fifth order are formed by
[`reparam_map_order5()`](https://statmodels7.github.io/distributions7/reference/reparam_map_order5.md)
from the fourth-order tables, which `map_derivs` supplies in closed
form.
