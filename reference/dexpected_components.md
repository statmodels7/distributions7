# The Components of a Derivative of the Expected Information

Enumerates, for a derivative order, the pairs \\(a,b)\\ and the
multisets of differentiating parameters, with the key under which
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and its relatives name each component.

## Usage

``` r
dexpected_components(params, m)
```

## Arguments

- params:

  The parameter names.

- m:

  The derivative order, 0 to 4.

## Value

A list of components, each with `key`, `a`, `b` and `S`.

## See also

[`dexpected_by_method()`](https://statmodels7.github.io/distributions7/reference/dexpected_by_method.md)
