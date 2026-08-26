# The Location Vector a Parameter List Describes

Returns the location of a multivariate distribution as a numeric vector
of length \\p\\. The parameters of a multivariate distribution are
scalars, so that every generic of the package can index them; this
generic puts the location back into the shape a reader thinks in.

## Usage

``` r
mv_location(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters, each component a single number.
  Aligned by the generic before dispatch.

## Value

A numeric vector of length \\p\\, named `v1`, ..., `vp` after the
coordinates of the response.

## Details

Not every multivariate family HAS a location: a Dirichlet is described
by concentrations and a Wishart by a scale and a count. The base-class
method therefore signals an error. Handing back the first \\p\\
parameters under a name that does not fit them would be worse. The
elliptical families register
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md),
which does exactly that and is right for them.

The location is the center of symmetry of the density. It is the MEAN as
well for a gaussian, and for a Student t only above one degree of
freedom; [`base::mean()`](https://rdrr.io/r/base/mean.html) is the
generic that answers about the moment.

## Notation

\\\mu\\ is the location and \\p\\ the dimension of one observation.

## See also

[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the matrix,
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md)
for the implementation the elliptical families use,
[`base::mean()`](https://rdrr.io/r/base/mean.html) for the moment, and
[`mv_location.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.multivariate_distrib.md)
for the refusal.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
mv_location(d, theta)
#> v1 v2 
#>  1 -1 

# A Student t has a location at every nu, and a mean only above nu = 1.
t2 <- mvstudent_t_distrib(2)
th <- c(theta, list(nu = 0.8))
rbind(location = mv_location(t2, th), mean = mean(t2, th))
#>           v1  v2
#> location   1  -1
#> mean     NaN NaN

# A Dirichlet has neither, and says so.
try(mv_location(dirichlet_distrib(3),
                list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)))
#> [1] 0.4260125 0.2583897 0.3155978
```
