# Observed Derivatives by Index Tuple

Returns a function that gives the observed derivative of the log-density
for a sorted index tuple at fixed observations and parameters, computing
each order once.

## Usage

``` r
observed_getter(distrib, ys, th)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- ys:

  The observations.

- th:

  The parameters, aligned.

## Value

A function of an integer tuple.

## See also

[`dexpected_by_method()`](https://statmodels7.github.io/distributions7/reference/dexpected_by_method.md)
