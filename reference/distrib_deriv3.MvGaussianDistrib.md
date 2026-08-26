# Multivariate Gaussian Third Derivatives

Computes every third derivative of the log-density in the parameters, in
closed form, on the matrix parametrization's own `param_d3()`. Writing
\\r = y - \mu\\ and \\P_t\\ for the precision's derivative array over
the matrix indices of the tuple: three mean indices give exactly zero,
the quadratic form being quadratic in \\\mu\\; two give \\-P_t\[i,
j\]\\; one gives \\(r P_t)\_i\\; none gives
\\\pm\tfrac{1}{2}\\\partial^3\log\lvert M\rvert - \tfrac{1}{2} r^\top
P_t r\\. The arrays \\P_t\\ come from the parametrization directly under
a precision parametrization, and from the expansion of the derivative of
an inverse under a covariance one.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. With
  `expected = TRUE` only its row count is used.

- theta:

  A named list of parameters, each component a single number.

- expected:

  Logical of length 1. When `TRUE` the expectation of each component is
  returned, by sampling. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. Every link of this family is the identity, so the two
  coincide.

- approx:

  Ignored: sampling is the only multivariate route to an expectation
  here. Present so that the signature matches the generic's.

- nsim:

  The number of draws used when `expected = TRUE`. Defaults to `10000`.
  Ignored otherwise.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 3)`. With `expected = TRUE` every vector is
constant.

## Details

With `expected = TRUE` the expectation is taken by SAMPLING: `nsim`
draws are made from the family and the observed components averaged over
them. There is no exact route here and no quadrature, so `approx` is not
read and the result carries Monte Carlo error of order `nsim^(-1/2)`.
Set a seed before calling if the result must be reproducible. A
component that does not depend on the response is returned exactly,
whatever `nsim` is, because averaging a constant returns it.

## Notation

\\\mu\\ is the mean, \\M\\ the matrix the parametrization carries,
\\\eta\\ its free vector, \\r = y - \mu\\ the centered response and
\\P_t\\ the precision's derivative array over the multiset \\t\\.

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
d <- mvgaussian_distrib(2)
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

# Against one stencil on the analytic Hessian.
h <- 1e-4
tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
c(exact = sum(d3[["mu1_mu2_sigma_L2.1"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
             sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
#>    exact  stencil 
#> 5.399435 5.399435 

# The expected version samples, so a component that carries the response
# moves with the seed while one that does not is exact.
set.seed(9)
e3 <- distrib_deriv3(d, y, theta, expected = TRUE, nsim = 4000)
c(sampled = e3[["mu1_mu1_sigma_L2.1"]][1],
  observed = mean(d3[["mu1_mu1_sigma_L2.1"]]))
#>    sampled   observed 
#> -0.9771222 -0.9771222 
```
