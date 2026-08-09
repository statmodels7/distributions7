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

For a subset \\A\\ of the coordinates the marginal density is

\$\$f_A(y_A; \theta) = \int\_{\mathbb{R}^{p - \lvert A \rvert}} f(y;
\theta)\\ \mathrm{d}y\_{A^{c}},\$\$

which for a Gaussian is Gaussian with \\\mu_A\\ and \\\Sigma\_{AA}\\,
and for a Student t is a Student t with \\\mu_A\\, \\\Sigma\_{AA}\\ and
the **same** \\\nu\\ (conditioning on the mixing variable of the
scale-mixture representation leaves a Gaussian, whose marginal is
Gaussian, and the mixture is then taken back). For a Dirichlet the
marginal of one coordinate is Beta with the same concentration \\\phi\\,
and for a multinomial it is binomial.

A marginal is not available in general: integrating a density over the
coordinates one is not interested in has no closed form for most
families. Every family the package ships has one, the elliptical ones by
subsetting the mean and the matrix and the two simplex families by the
identities above. For a family without one the generic signals an error
rather than approximating, since a quadrature over the discarded
coordinates would be a different object under the same name.

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
