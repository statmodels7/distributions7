# A Multivariate t Proposal for Integrating the Density

Draws from a multivariate t with the family's location, twice its scale
matrix and half its degrees of freedom, and returns the log-density of
the draws under that law. Its tails are heavier than the family's at
every \\\nu\\, so the importance ratio is bounded and its second moment
finite; the default gaussian proposal needs a finite covariance, which a
t does not have at \\\nu \le 2\\.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- theta:

  A named list of parameters, already aligned by the generic.

- n:

  The number of draws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `y`, an \\n \times p\\ numeric matrix, and `logd`, a
numeric vector of length \\n\\.

## See also

[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
for the generic and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
for the consumer.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5, nu = 1.5)
set.seed(1)
r <- mv_reference_draw(d, theta, 20000)
mean(exp(distrib_pdf(d, r$y, theta, log = TRUE) - r$logd))
#> [1] 0.9949695
```
