# Require Scalar Parameters, and Flatten Them

Rejects a `theta` whose components are not single numbers, and returns
the rest as one named numeric vector in `distrib@params` order. The
families here take ONE parameter value for the whole sample: a
parametrization describes one matrix, not one per row, and a parameter
that varies by observation belongs to a model, which is the layer above
this one.

## Usage

``` r
mv_flat_theta(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned by
  [`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md).
  Every component must have length 1.

## Value

A named numeric vector of length `distrib@n_params`, named and ordered
as `distrib@params`.

## Details

A univariate family reads a `theta` component of length \\n\\ as one
value per observation, so without this check a caller who wrote that
here would get silent recycling against the \\p\\ columns instead of an
error. The message names the offending parameter.

## See also

[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
which normalizes the names first, and
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
for the convention.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
distributions7:::mv_flat_theta(d, distributions7:::align_theta(d, theta))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1 
#>          0.0          0.0          0.0          0.0          0.5 

# A component longer than one is rejected by name, where a univariate
# family would read it as one value per observation.
bad <- theta
bad$mu1 <- c(0, 1)
try(distributions7:::mv_flat_theta(d, bad))
#> Error : Parameter 'mu1' has length 2. A multivariate distribution takes one
#>   value per parameter for the whole sample; parameters that vary by
#>   observation belong to a model.
```
