# A Marginal of a Multivariate Distribution

Returns the distribution of a subset of the coordinates, together with
its parameters. For a subset \\A\\, \$\$f_A(y_A; \theta) =
\int\_{\mathbb{R}^{\lvert A^{c}\rvert}} f(y; \theta)\\
\mathrm{d}y\_{A^{c}},\$\$ which the shipped families have in closed
form: any subset for the gaussian and the Student t, and a single
coordinate for the Dirichlet and the multinomial. The method of the
class signals an error, because a quadrature over the discarded
coordinates would give a different object under the same name.

## Usage

``` r
mv_marginal(distrib, theta, which, ...)
```

## Arguments

- distrib:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters, each component a single number.
  Aligned by the generic before dispatch.

- which:

  A vector of distinct whole numbers in `1:distrib@n_dim`, the
  coordinates to keep. The generic checks it before dispatch and signals
  an error naming the range for an empty, missing, non-integer,
  duplicated or out-of-range value. The Dirichlet and the multinomial
  accept a single coordinate.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A named list with `distrib`, the marginal distribution object, and
`theta`, its parameters as a named list.

## What each family returns

A gaussian's marginal is a gaussian with \\\mu_A\\ and \\\Sigma\_{AA}\\.
A Student t's is a Student t with \\\mu_A\\, \\\Sigma\_{AA}\\ and the
same \\\nu\\: conditionally on the mixing variable of the scale-mixture
representation the vector is gaussian, whose marginal is gaussian, and
mixing again gives a Student t. A Dirichlet's marginal of one coordinate
is a beta with the same concentration \\\phi\\, and a multinomial's is a
binomial. For several coordinates these two signal an error, because the
subvector has the family's form only once a coordinate for the remaining
mass is added.

## The returned object

A new distribution of the reduced dimension, whose parameters are its
own. The elliptical families return an object of their own class on an
unstructured log-Cholesky matrix, carrying the covariance or the
precision as the parent does, so a gaussian's one-coordinate marginal is
a one-dimensional `MvGaussian1Distrib` and signals an error for the
distribution function. The simplex-valued families return univariate
objects, a beta and a binomial, which have every method of a univariate
family.

## The plot

A panel of the pairs plot of
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)
is a marginal of one or two coordinates, so the plot is available for
the gaussian and the Student t, and signals the error above for the
Dirichlet and the multinomial.

## Notation

\\A\\ is the retained subset of coordinates, \\A^c\\ its complement,
\\\mu\\ the location, \\\Sigma\\ the matrix that the family carries,
\\\nu\\ a Student t's degrees of freedom and \\\phi\\ a Dirichlet's
concentration.

## See also

[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
whose panels are these marginals, and
[`mv_marginal.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvGaussianDistrib.md),
[`mv_marginal.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvStudentTDistrib.md)
and
[`mv_marginal.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.multivariate_distrib.md)
for the methods and the error.

## Examples

``` r
d <- mvgaussian1_distrib(3)
theta <- as.list(stats::setNames(
  c(1, 2, 3, 0, 0, 0, 0.3, 0.2, 0.1), d@params))

# The marginal of two coordinates is a two-dimensional gaussian on the
# corresponding block of the covariance.
m <- mv_marginal(d, theta, c(1, 2))
mv_sigma(m$distrib, m$theta)
#>     v1   v2
#> v1 1.0 0.30
#> v2 0.3 1.09
mv_sigma(d, theta)[1:2, 1:2]
#>     v1   v2
#> v1 1.0 0.30
#> v2 0.3 1.09

# A Student t's marginal keeps the same degrees of freedom.
t3 <- mvstudent_t1_distrib(3)
th <- as.list(stats::setNames(c(unlist(theta), 5), t3@params))
c(full = th$nu, marginal = mv_marginal(t3, th, c(1, 3))$theta$nu)
#>     full marginal 
#>        5        5 

# A Dirichlet's is a beta with the same concentration, a univariate object.
b <- mv_marginal(dirichlet_distrib(3),
                 list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
class(b$distrib)[1]
#> [1] "distributions7::Beta1Distrib"
b$theta
#> $mu
#> [1] 0.4260125
#> 
#> $phi
#> [1] 8
#> 

# 'which' is checked before dispatch.
try(mv_marginal(d, theta, c(1, 1)))
#> Error : 'which' must be distinct coordinates in 1:3.
try(mv_marginal(d, theta, 4))
#> Error : 'which' must be distinct coordinates in 1:3.
```
