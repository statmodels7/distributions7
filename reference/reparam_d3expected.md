# The Third Derivative of the Expected Information of a Reparametrized Family

Reads the parent's expected information and its first three derivatives
at the mapped parameters and the map's tables to the fourth order, and
passes them to
[`dexpected_chain3()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain3.md).

## Usage

``` r
reparam_d3expected(distrib, y, theta, threads = 1L)
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
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md).
