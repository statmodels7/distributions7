# Scale Standard Deviations and Correlations of a Multivariate t

Returns the square roots of the diagonal of the SCALE matrix and the
correlations it implies, with the closed-form Jacobian
[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
supplies. The correlations are the response's as well: the covariance is
\\\nu\Sigma/(\nu-2)\\ where it exists, and a positive multiple of a
matrix leaves its correlations alone. The diagonal quantities are NOT
standard deviations of the response, and are named `scale_sd_v1`, ...,
`scale_sd_vp` and blocked as `"Scale standard deviations"` to say so.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- theta:

  A named list of parameters, already aligned by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `value`, `jacobian`, `transform` and `block`, as
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
documents: \\p\\ scale standard deviations on the log scale and
\\p(p-1)/2\\ correlations on Fisher's \\z\\, plus whatever the matrix
parametrization declares.

## Details

The degrees of freedom appear in no quantity here. \\\nu\\ is already an
interpretable parameter on its own scale, so
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
reports it and this generic does not repeat it. That is also why the
\\\nu\\ column of the Jacobian is zero.

## Notation

\\\Sigma\\ is the scale matrix, \\\nu\\ the degrees of freedom and \\p\\
the dimension.

## See also

[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
for the closed-form Jacobian,
[`mv_derived.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md),
whose diagonal quantities are genuine standard deviations,
[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
for the covariance, and
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5, nu = 6)
der <- mv_derived(d, theta)
der$value
#> scale_sd_v1 scale_sd_v2   cor_v1_v2 
#>   1.0000000   1.1180340   0.4472136 
der$block
#>                 scale_sd_v1                 scale_sd_v2 
#> "Scale standard deviations" "Scale standard deviations" 
#>                   cor_v1_v2 
#>              "Correlations" 

# The diagonal quantities are the scale matrix's, not the response's: the
# response is more spread than they say by sqrt(nu / (nu - 2)).
c(scale_sd = der$value[["scale_sd_v1"]],
  response_sd = sqrt(variance(d, theta)[1, 1]))
#>    scale_sd response_sd 
#>    1.000000    1.224745 

# The correlation is the same either way, the two matrices being positive
# multiples of each other.
cv <- function(m) m[1, 2] / sqrt(m[1, 1] * m[2, 2])
c(reported = der$value[["cor_v1_v2"]],
  from_covariance = cv(variance(d, theta)))
#>        reported from_covariance 
#>       0.4472136       0.4472136 

# nu enters no quantity here, so its Jacobian column is zero.
der$jacobian[, "nu"]
#> scale_sd_v1 scale_sd_v2   cor_v1_v2 
#>           0           0           0 
```
