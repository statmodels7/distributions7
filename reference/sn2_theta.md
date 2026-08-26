# The Direct Parameters a Centered Triple Implies

Takes the sign of the skewness off its plain value and runs
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md),
returning the location, scale and shape that
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
takes. Every probability function of the centered family calls this and
then delegates.

## Usage

``` r
sn2_theta(theta)
```

## Arguments

- theta:

  A list with `mu`, `sigma` and `gamma1`, in that order, each a numeric
  vector. It is read positionally, so it must already be aligned; the
  methods that call it have been through
  [`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md).

## Value

A named list with `mu`, `sigma` and `alpha`: the location \\\xi\\, the
scale \\\omega\\ and the shape \\\alpha\\. The names are the parent's,
so the result passes straight into
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)'s
methods.

## See also

[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
for the map itself and
[`sn2_chain()`](https://statmodels7.github.io/distributions7/reference/sn2_chain.md)
for the derivative route that uses the same map.

## Examples

``` r
distributions7:::sn2_theta(list(mu = 0, sigma = 1, gamma1 = 0.5))
#> $mu
#> [1] -1.052209
#> 
#> $sigma
#> [1] 1.451601
#> 
#> $alpha
#> [1] 2.173758
#> 

# A negative skewness reflects the shape and moves the location the other way.
distributions7:::sn2_theta(list(mu = 0, sigma = 1, gamma1 = -0.5))
#> $mu
#> [1] 1.052209
#> 
#> $sigma
#> [1] 1.451601
#> 
#> $alpha
#> [1] -2.173758
#> 
```
