# A Few Points of the Support

Three quantiles of the family at the probe parameters, falling back on
the interior of its bounds where the quantile function refuses.

## Usage

``` r
.kink_probe_y(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameter values.

## Value

A numeric vector.

## Examples

``` r
d <- laplace_distrib()
distributions7:::.kink_probe_y(d, distributions7:::.kink_probe_theta(d))
#> [1] -0.6931472  0.0000000  0.6931472
```
