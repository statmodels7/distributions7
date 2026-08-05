# Generate Random Parameters

Generates sensible random parameters for a distribution based on its
mathematical domain.

## Usage

``` r
generate_random_theta(distrib, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of parameter values, one per element of `distrib@params`,
each inside that parameter's bounds.

## Examples

``` r
set.seed(1)
generate_random_theta(gamma2_distrib())
#> $mu
#> [1] 1.666501
#> 
#> $sigma2
#> [1] 2.295531
#> 
```
