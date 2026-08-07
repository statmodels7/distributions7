# The Mean Vector and Covariance a Parameter List Describes

Assembles the mean vector and the covariance matrix of a multivariate
distribution from its flat parameter list.

## Usage

``` r
mv_location(distrib, theta)

mv_sigma(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list or vector of parameters.

  Both are generics whose base-class method rejects: not every
  multivariate family has a location, and one that does not should say
  so rather than hand back its first p parameters under a name that does
  not fit them.

## Value

A numeric vector of length \\p\\ for `mv_location()`, and a \\p \times
p\\ matrix for `mv_sigma()`.

## Details

The parameters of a multivariate distribution are scalars, so that every
generic of the package can index them, and these two functions put them
back into the shapes a reader thinks in. `mv_sigma()` returns the matrix
the PARAMETRIZATION carries, whichever side the matrix parameter
describes: the covariance for a gaussian, and the scale matrix for a
Student t, whose covariance is \\\nu\Sigma/(\nu-2)\\ and does not exist
below two degrees of freedom. The moment is
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md),
and keeping the two apart is what lets a heavy-tailed family be
described at all.

## See also

[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)
mv_location(d, theta)
#> v1 v2 
#>  1 -1 
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25
```
