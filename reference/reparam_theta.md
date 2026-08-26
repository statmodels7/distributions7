# The Parent's Parameters at the New Ones

Runs the map on plain numbers and returns the parent's parameters in the
parent's own order, which is the form every probability function needs
before delegating. The new parameters are aligned and validated first,
so a missing or out-of-bounds component throws here and not several
frames down.

## Usage

``` r
reparam_theta(distrib, theta)
```

## Arguments

- distrib:

  A reparametrized distribution.

- theta:

  A named list of the new parameters, on the new parameter scale.
  Components may be vectors.

## Value

A named list of the parent's parameters, in the parent's order.

## Details

A map that fails to return some parent parameter is caught with a
message naming the parameter and the parent's full list. The check is
cheap and is repeated on every call: a map is ordinary R and may return
different names at different points, so the construction cannot settle
it once.

## See also

[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
for the map's derivatives;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md);
[`ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md).
