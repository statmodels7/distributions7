# Matrix Entries as the Default Interpretable Quantities

Returns the distinct entries of the matrix that
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
produces, named `sigma_vi_vj` after their coordinates, with a Jacobian
from one central difference in each parameter, followed by the block
that
[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md)
returns when the parametrization declares one. It is the method of a
family that registers none of its own, among them the Dirichlet and the
multinomial, whose matrix is the covariance of the response. The entries
are variances and covariances, where a closed-form method reports
standard deviations and correlations.

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
documents. There are `p * (p + 1) / 2` matrix entries, blocked as
`"Variances"` and `"Covariances"`, followed by the quantities that the
matrix parametrization declares through `param_readable()` when it
declares any.

## Details

The diagonal entries are variances and are intervalled on the log scale;
the off-diagonal ones are unconstrained given the diagonal and are
intervalled on their own scale. The closed-form methods report
correlations on Fisher's \\z\\, whose intervals stay inside \\(-1, 1)\\;
the interval of a covariance has no such bound.

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
d <- mvgaussian1_distrib(2)
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

# The entries are the lower triangle of the matrix.
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25

# Against the closed-form method, which reports standard deviations and a
# correlation in place of variances and a covariance.
mv_derived(d, theta)$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 

# The Dirichlet uses this method: its covariance, then its probabilities.
dd <- dirichlet_distrib(3)
mv_derived(dd, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8))$block
#>     sigma_v1_v1     sigma_v2_v1     sigma_v3_v1     sigma_v2_v2     sigma_v3_v2 
#>     "Variances"   "Covariances"   "Covariances"     "Variances"   "Covariances" 
#>     sigma_v3_v3              p1              p2              p3 
#>     "Variances" "Probabilities" "Probabilities" "Probabilities" 
```
