# Interpretable Quantities of a Multivariate Distribution

Returns the quantities a reader of a multivariate fit actually wants —
standard deviations, correlations, and whatever else the family's matrix
parameter means — together with the Jacobian needed to carry standard
errors onto them.

## Usage

``` r
mv_derived(distrib, theta, ...)
```

## Arguments

- distrib:

  An object inheriting from class
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters.

- ...:

  Passed to methods.

## Value

A list with

- `value`:

  a named numeric vector of the quantities;

- `jacobian`:

  a matrix with one row per quantity and one column per parameter of
  `distrib`;

- `transform`:

  a character vector, one of `"identity"`, `"log"` or `"atanh"` per
  quantity, naming the scale its interval is built on;

- `block`:

  a character vector labeling the group each quantity belongs to, used
  to lay the printed summary out.

## Details

The free values of a parameters7 structure are coordinates chosen so
that an optimizer can move freely; they are not quantities anyone reads.
The logarithm of the third diagonal entry of a Cholesky factor has an
estimate and a standard error, and neither answers a question. This
generic names the quantities that do, and supplies \\\partial
g/\partial\theta\\ so that
[`mv_summary`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
can apply the delta method to them.

Each quantity also declares the scale its confidence interval should be
built on, exactly as
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
builds a univariate interval on the link scale and maps it back: a
standard deviation is intervalled on the log scale so that the interval
cannot reach zero, a correlation on Fisher's \\z =
\mathrm{artanh}(\rho)\\ so that it cannot leave \\(-1, 1)\\, and an
unconstrained quantity on its own scale.

The default method, registered on
[`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md),
returns the distinct entries of the matrix
[`mv_sigma`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
produces, named after the coordinates they belong to, with a Jacobian
obtained by one central difference. A family whose matrix is not a
covariance therefore still reports something on its original scale,
which is better than reporting a Cholesky coordinate.

## See also

[`mv_summary`](https://statmodels7.github.io/distributions7/reference/mv_summary.md),
[`mv_sigma`](https://statmodels7.github.io/distributions7/reference/mv_location.md)

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
mv_derived(d, theta)$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 
```
