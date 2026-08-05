# Fix Parameters of a Distribution at Known Values

Returns the distribution obtained by holding some parameters of
`distrib` at known values, leaving only the others to be supplied and
estimated.

## Usage

``` r
fixed(distrib, ...)
```

## Arguments

- distrib:

  The distribution whose parameters are to be fixed.

- ...:

  The fixed values, named after the parameters they fix, as in
  `fixed(gaussian1_distrib(), mu = 0)`.

## Value

An object of class `FixedContinuousDistrib` or `FixedDiscreteDistrib`,
matching the parent.

## Details

The result is the same law with a smaller parameter set: `theta` carries
only the free parameters, every generic answers as the parent does at
the full vector, and the derivative components are the parent's
restricted to the free indices – a subvector of the score, a submatrix
of the Hessian, sub-arrays at orders three and four. Nothing is
recomputed and no normalizing constant changes, so the parent's closed
forms are used throughout, and
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
estimates the free parameters with standard errors and intervals for
them alone.

Fixed values are single numbers, strictly inside the open domain of
their parameter. Fixing a parameter of a distribution that is already a
fixed-parameter wrapper collapses the two into one wrapper around the
original parent. Fixing every parameter is allowed and gives a fully
known distribution with an empty parameter set.

The per-parameter smoothness declaration travels with the free
parameters, so fixing the location of a Laplace distribution leaves a
distribution whose remaining parameter is smooth.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
# a gaussian with known mean: only sigma remains
d <- fixed(gaussian1_distrib(), mu = 0)
d@params
#> [1] "sigma"

theta <- list(sigma = 2)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1760327 0.1994711 0.1760327
distrib_gradient(d, c(-1, 0, 1), theta)
#> $sigma
#> [1] -0.375 -0.500 -0.375
#> 

# the score is the corresponding component of the parent's
full <- distrib_gradient(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 2))
all.equal(distrib_gradient(d, c(-1, 0, 1), theta)$sigma, full$sigma)
#> [1] TRUE

# fixing everything gives a fully known distribution
d0 <- fixed(gaussian1_distrib(), mu = 0, sigma = 1)
distrib_pdf(d0, 0, list())
#> [1] 0.3989423
```
