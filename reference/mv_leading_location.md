# The First p Parameters, Read as a Location

Reads the first \\p\\ components of `theta` as the location vector. It
is the implementation with which both elliptical families register
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md).
It is correct where the leading \\p\\ parameters of a family are its
location, which is a property of the family and not of the shape of
`theta`.

## Usage

``` r
mv_leading_location(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object whose leading \\p\\ parameters are the location.

- theta:

  A named list of parameters, already aligned.

## Value

A numeric vector of length `distrib@n_dim`, named `v1`, ..., `vp`.

## See also

[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic and
[`mv_location.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md)
for a registration that uses this.

## Examples

``` r
d <- mvgaussian1_distrib(3)
theta <- as.list(stats::setNames(c(1, -2, 0.5, rep(0, 6)), d@params))
th <- distributions7:::align_theta(d, theta)
distributions7:::mv_leading_location(d, th)
#>   v1   v2   v3 
#>  1.0 -2.0  0.5 

# It is what mv_location() returns for this family, the method being this
# function.
mv_location(d, theta)
#>   v1   v2   v3 
#>  1.0 -2.0  0.5 
```
