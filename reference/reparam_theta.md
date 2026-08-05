# The Parent's Parameters at the New Ones

Runs the map on plain numbers, which is what every probability function
needs before delegating to the parent.

## Usage

``` r
reparam_theta(distrib, theta)
```

## Arguments

- distrib:

  A reparametrized distribution.

- theta:

  A named list of the new parameters.

## Value

A named list of the parent's parameters.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
