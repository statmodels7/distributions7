# Location of a Multivariate Student t

Returns the location \\\mu\\, the first \\p\\ parameters read off
`theta` in order. The method is
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md),
shared with the gaussian. For this family the location is the center of
symmetry at every admissible \\\nu\\ and the mean only for \\\nu \> 1\\,
which is why the generic is named for a location:
[`mean.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MvStudentTDistrib.md)
returns `NaN` below one degree of freedom and this returns the location
regardless.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- theta:

  A named list of parameters, one number each. Only the \\p\\ location
  components are read.

## Value

A numeric vector of length \\p\\, named `v1`, ..., `vp` after the
coordinates of the response.

## See also

[`mean.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MvStudentTDistrib.md)
for the moment,
[`mv_sigma.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
for the matrix, and
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

mv_location(d, theta)
#>   v1   v2 
#>  0.5 -0.3 

# Above one degree of freedom the location is the mean; below it, only the
# location survives.
t2 <- theta; t2$nu <- 0.8
rbind(location = mv_location(d, t2), mean = mean(d, t2))
#>           v1   v2
#> location 0.5 -0.3
#> mean     NaN  NaN
```
