# A Trial Parameter Value Inside Every Interval

The midpoint of a bounded interval, one unit inside a half-line, and
zero on the whole line. The same rule
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
probes a map with.

## Usage

``` r
.kink_probe_theta(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

A named list, one entry per parameter.

## Examples

``` r
distributions7:::.kink_probe_theta(laplace_distrib())
#> $mu
#> [1] 0
#> 
#> $sigma
#> [1] 1
#> 
```
