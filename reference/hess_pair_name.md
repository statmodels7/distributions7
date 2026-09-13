# The Hessian Name of a Pair of Parameters

`params[a]_params[b]` in whichever order
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
lists it.

## Usage

``` r
hess_pair_name(params, a, b)
```

## Arguments

- params:

  A character vector of parameter names.

- a, b:

  Indices into `params`.

## Value

A single string.
