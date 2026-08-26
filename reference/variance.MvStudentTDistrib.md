# Covariance of a Multivariate Student t

Returns \\\operatorname{Var}(Y) = \nu\Sigma/(\nu-2)\\ for \\\nu \> 2\\,
and a matrix of `Inf` otherwise. The scale matrix \\\Sigma\\ is what the
parametrization carries and exists at every admissible \\\nu\\; the
covariance is a moment, and the two differ by that factor, which is 3 at
\\\nu = 3\\ and approaches 1 as the tail lightens.
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
returns the first and this returns the second.

## Arguments

- x:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- theta:

  A named list of parameters, each component a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A \\p \times p\\ numeric matrix with both dimnames `v1`, ..., `vp`: the
covariance for \\\nu \> 2\\, and `Inf` in every entry for \\\nu \le 2\\,
the boundary \\\nu = 2\\ included.

## See also

[`mv_sigma.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
for the scale matrix,
[`mean.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MvStudentTDistrib.md)
for the first moment, which needs only \\\nu \> 1\\, and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

variance(d, theta)
#>           v1        v2
#> v1 1.8321041 0.6631026
#> v2 0.6631026 1.2454801
mv_sigma(d, theta)
#>           v1        v2
#> v1 1.2214028 0.4420684
#> v2 0.4420684 0.8303200
all.equal(variance(d, theta), (6 / 4) * mv_sigma(d, theta))
#> [1] TRUE

# The factor blows up as nu falls to 2 and the covariance goes with it.
vapply(c(10, 3, 2.5, 2.05, 2), function(nu) {
  t2 <- theta; t2$nu <- nu
  variance(d, t2)[1, 1]
}, numeric(1))
#> [1]  1.526753  3.664208  6.107014 50.077513       Inf

# Correlations are unaffected: a positive multiple leaves them alone.
cv <- function(m) m[1, 2] / sqrt(m[1, 1] * m[2, 2])
c(scale = cv(mv_sigma(d, theta)), covariance = cv(variance(d, theta)))
#>      scale covariance 
#>  0.4389724  0.4389724 
```
