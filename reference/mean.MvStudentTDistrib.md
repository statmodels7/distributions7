# Mean of a Multivariate Student t

Returns \\\mathbb{E}\[Y\] = \mu\\ for \\\nu \> 1\\, and a vector of
`NaN` otherwise. For one degree of freedom or fewer the defining
integral does not converge: the density decays like \\\lVert
y\rVert^{-(\nu+p)}\\ and the first absolute moment integrates \\\lVert
y\rVert^{1-\nu-p}\\ over a shell of surface \\\lVert y\rVert^{p-1}\\,
which is finite only for \\\nu \> 1\\. The location \\\mu\\ is a
parameter and exists at every \\\nu\\, while the mean is a moment; for
\\\nu \le 1\\ this method returns `NaN` and
[`mv_location.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvStudentTDistrib.md)
returns \\\mu\\.

## Arguments

- x:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- theta:

  A named list of parameters; each component is a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length \\p\\, named `v1`, ..., `vp`, holding the
location for \\\nu \> 1\\ and `NaN` throughout for \\\nu \le 1\\.

## See also

[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
for the second moment, which needs \\\nu \> 2\\,
[`mv_location.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvStudentTDistrib.md)
for the location, which exists at every \\\nu\\, and
[`base::mean()`](https://rdrr.io/r/base/mean.html) for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

mean(d, theta)
#>   v1   v2 
#>  0.5 -0.3 

# Below one degree of freedom the mean does not exist, while the location
# and the density do.
t2 <- theta; t2$nu <- 0.8
mean(d, t2)
#>  v1  v2 
#> NaN NaN 
mv_location(d, t2)
#>   v1   v2 
#>  0.5 -0.3 
distrib_pdf(d, c(0, 0), t2)
#> [1] 0.08456772
```
