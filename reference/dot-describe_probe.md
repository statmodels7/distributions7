# Describe a Trial Parameter Value

Renders the point a map was probed at as `(a = 1, b = 2)`, so that a
failure in
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)'s
construction check names the point instead of leaving the caller to
guess it.

## Usage

``` r
.describe_probe(probe)
```

## Arguments

- probe:

  A named list of parameter values, each of length 1.

## Value

A single string.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
the one caller.
