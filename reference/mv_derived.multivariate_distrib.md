# Matrix Entries as the Default Interpretable Quantities

Returns the distinct entries of the matrix
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
produces, named `sigma_vi_vj` after the coordinates they belong to, with
a Jacobian from one central difference in each parameter. This is the
fallback for a family that registers no method of its own: reporting the
matrix on its own scale is worth more to a reader than reporting a
Cholesky coordinate, even where the entries are not the standard
deviations and correlations a closed-form method would give.

## Arguments

- distrib:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  that has no method of its own.

- theta:

  A named list of parameters, already aligned by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `value`, `jacobian`, `transform` and `block`, as
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
documents. There are `p * (p + 1) / 2` quantities, blocked as
`"Variances"` and `"Covariances"`.

## Details

The diagonal entries are variances and are intervalled on the log scale;
the off-diagonal ones are unconstrained given the diagonal and are
intervalled on their own. That is weaker than the closed-form methods,
whose correlations ride Fisher's \\z\\ and so cannot leave \\(-1, 1)\\.

## See also

[`mv_derived.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md)
and
[`mv_derived.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvStudentTDistrib.md)
for the two closed-form methods,
[`mv_entry_index()`](https://statmodels7.github.io/distributions7/reference/mv_entry_index.md)
for the labeling, and
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)

# The gaussian overrides, so reach the default method directly.
base <- S7::method(mv_derived, multivariate_distrib)
d0 <- base(d, distributions7:::align_theta(d, theta))
d0$value
#> sigma_v1_v1 sigma_v2_v1 sigma_v2_v2 
#>        1.00        0.50        1.25 
d0$transform
#> sigma_v1_v1 sigma_v2_v1 sigma_v2_v2 
#>       "log"  "identity"       "log" 

# The entries really are the matrix, read off its lower triangle.
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25

# Against the closed-form method, which reports standard deviations and a
# correlation instead of variances and a covariance.
mv_derived(d, theta)$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 
```
