# Multivariate Gaussian Third Derivatives

Computes every third derivative of the log-density in the parameters, in
closed form. Writing \\r = y - \mu\\ and \\P_t\\ for the derivative of
the precision over the matrix indices of the tuple, a tuple with three
mean indices gives exactly zero, because the log-density is quadratic in
\\\mu\\; a tuple with two gives \\-P_t\[i, j\]\\; a tuple with one gives
\\(r P_t)\_i\\; and a tuple without mean indices gives
\\\pm\tfrac{1}{2}\\\partial^3\log\lvert M\rvert - \tfrac{1}{2} r^\top
P_t r\\. The arrays \\P_t\\ come from
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md).

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. With
  `expected = TRUE` only its row count is used.

- theta:

  A named list of parameters in which each component is a single number.

- expected:

  Logical of length 1. When `TRUE` the expectation of each component is
  returned. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.
  Every link of this family is the identity, so the two coincide.

- approx:

  Not read: the expectation is exact.

- nsim:

  Not read, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 3)`. With `expected = TRUE` every vector is
constant.

## Details

With `expected = TRUE` the expectation is exact: the response enters
through \\r\\, whose expectation is zero, and through \\r^\top P_t r\\,
whose expectation is \\\operatorname{tr}(P_t\Sigma)\\. `approx` and
`nsim` are not read.

## Notation

\\\mu\\ is the mean, \\M\\ the matrix that the parametrization carries,
\\\Sigma\\ the covariance, \\\eta\\ the free vector, \\r = y - \mu\\ the
centered response and \\P_t\\ the derivative of the precision over the
multiset \\t\\.

## See also

[`distrib_deriv4.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MvGaussianDistrib.md)
for the next order,
[`distrib_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
for the second,
[`mvg_higher()`](https://statmodels7.github.io/distributions7/reference/mvg_higher.md)
for the shared engine, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

d3 <- distrib_deriv3(d, y, theta)
length(d3)
#> [1] 35

# Any component naming three means is exactly zero.
d3[["mu1_mu1_mu2"]]
#> [1] 0 0 0 0

# The result is compared with one stencil on the analytic Hessian.
h <- 1e-4
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
c(exact = sum(d3[["mu1_mu2_sigma_L2.1"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
             sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
#>    exact  stencil 
#> 5.399435 5.399435 

# The expectation in closed form is compared with an average over a large
# sample.
e3 <- distrib_deriv3(d, y, theta, expected = TRUE)
set.seed(9)
big <- distrib_deriv3(d, distrib_rng(d, 50000, theta), theta)
c(closed = e3[["sigma_log_L1_sigma_log_L1_sigma_L2.1"]][1],
  sampled = mean(big[["sigma_log_L1_sigma_log_L1_sigma_L2.1"]]))
#>    closed   sampled 
#> -1.790190 -1.801327 
```
