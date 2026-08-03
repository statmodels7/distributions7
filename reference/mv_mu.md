# The Mean Vector and Covariance a Parameter List Describes

Assembles the mean vector and the covariance matrix of a multivariate
distribution from its flat parameter list.

## Usage

``` r
mv_mu(distrib, theta)

mv_sigma(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list or vector of parameters.

## Value

A numeric vector of length \\p\\ for `mv_mu()`, and a \\p \times p\\
matrix for `mv_sigma()`.

## Details

The parameters of a multivariate distribution are scalars, so that every
generic of the package can index them, and these two functions put them
back into the shapes a reader thinks in. `mv_sigma()` returns the
covariance whichever side the structure parametrises.

## See also

[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 1, mu2 = -1, log_L1 = 0, log_L2 = 0, L2.1 = 0.5)
mv_mu(d, theta)
#> v1 v2 
#>  1 -1 
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25
```
