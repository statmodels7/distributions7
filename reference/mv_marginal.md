# A Marginal of a Multivariate Distribution

Returns the distribution of a subset of the coordinates, together with
the parameters that describe it.

## Usage

``` r
mv_marginal(distrib, theta, which, ...)
```

## Arguments

- distrib:

  An object inheriting from class
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters.

- which:

  An integer vector of coordinates to keep.

- ...:

  Passed to methods.

## Value

A list with `distrib`, the marginal distribution object, and `theta`,
its parameters.

## Details

A marginal is not available in general: integrating a density over the
coordinates one is not interested in has no closed form for most
families. It is available for the elliptical ones, where the marginal
belongs to the same family with the mean and the matrix subsetted, and
those are the ones this generic has methods for. For a family without
one the generic signals an error rather than approximating, since a
quadrature over the discarded coordinates would be a different object
under the same name.

This is what makes a picture of a multivariate distribution possible at
all: a panel of a pairs plot shows a marginal, so the plot exists
exactly when the marginals do.

## See also

[`plot.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)

## Examples

``` r
d <- mvgaussian_distrib(3)
theta <- as.list(stats::setNames(c(1, 2, 3, 0, 0, 0, 0.3, 0.2, 0.1), d@params))

# the marginal of the first two coordinates is a two-dimensional gaussian
m <- mv_marginal(d, theta, c(1, 2))
mv_sigma(m$distrib, m$theta)
#>     v1   v2
#> v1 1.0 0.30
#> v2 0.3 1.09

# and it is the corresponding block of the full covariance
mv_sigma(d, theta)[1:2, 1:2]
#>     v1   v2
#> v1 1.0 0.30
#> v2 0.3 1.09
```
