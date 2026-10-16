# The Location Vector That a Parameter List Describes

Returns the location of a multivariate distribution as a numeric vector
of length \\p\\. The parameters of a multivariate distribution are
scalars, so that every generic of the package can index them, and this
generic returns the location in the shape of the response.

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
  The generic aligns it before dispatch.

## Value

A numeric vector of length \\p\\. The elliptical families name it `v1`,
..., `vp` after the coordinates of the response; the Dirichlet and the
multinomial return it without names.

## Details

The elliptical families return the center of symmetry of the density,
which is the mean for a gaussian and, for a Student t, the mean only
when \\\nu \> 1\\; they register
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md).
The Dirichlet returns its mean vector, and the multinomial its mean
vector of counts, the probabilities times the number of trials. The
method of the class signals an error for a family without a method of
its own, the first \\p\\ parameters of a family not being a location in
general. [`base::mean()`](https://rdrr.io/r/base/mean.html) is the
generic for the moment.

## Notation

\\\mu\\ is the location and \\p\\ the dimension of one observation.

## See also

[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the matrix,
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md)
for the implementation that the elliptical families use,
[`base::mean()`](https://rdrr.io/r/base/mean.html) for the moment, and
[`mv_location.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.multivariate_distrib.md)
for the method of the class.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
mv_location(d, theta)
#> v1 v2 
#>  1 -1 

# A Student t has a location at every nu, and a mean only above nu = 1.
t2 <- mvstudent_t1_distrib(2)
th <- c(theta, list(nu = 0.8))
rbind(location = mv_location(t2, th), mean = mean(t2, th))
#>           v1  v2
#> location   1  -1
#> mean     NaN NaN

# A Dirichlet returns its mean vector, a point of the simplex.
mv_location(dirichlet_distrib(3),
            list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8))
#> [1] 0.4260125 0.2583897 0.3155978
```
