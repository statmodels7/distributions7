# A Marginal of a Multivariate Distribution

Returns the distribution of a subset of the coordinates, together with
its parameters. For a subset \\A\\, \$\$f_A(y_A; \theta) =
\int\_{\mathbb{R}^{\lvert A^{c}\rvert}} f(y; \theta)\\
\mathrm{d}y\_{A^{c}},\$\$ which every family this package ships has in
closed form. A family without one is refused: a quadrature over the
discarded coordinates would be a different object under the same name.

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

  An integer vector of coordinates to keep. The generic checks it before
  dispatch: it must be non-empty, free of `NA`, free of duplicates, and
  inside `1:distrib@n_dim`. Anything else is an error naming the range.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list with `distrib`, the marginal distribution object, and
`theta`, its parameters as a named list.

## What each family gives

A gaussian's marginal is a gaussian with \\\mu_A\\ and \\\Sigma\_{AA}\\.
A Student t's is a Student t with \\\mu_A\\, \\\Sigma\_{AA}\\ and THE
SAME \\\nu\\: conditioning on the mixing variable of the scale-mixture
representation leaves a gaussian, whose marginal is gaussian, and the
mixture is then taken back. A Dirichlet's marginal of one coordinate is
a beta with the same concentration \\\phi\\, and a multinomial's is a
binomial.

## What the returned object is

A FRESH distribution of the reduced dimension, whose parameters are its
own. The elliptical families return an object of their own class on an
unstructured matrix, so a gaussian's one-coordinate marginal is a
one-dimensional `MvGaussianDistrib` and still refuses a distribution
function; the simplex-valued families return genuinely univariate
objects, a beta and a binomial, which answer everything a univariate
family does.

## Why the plot depends on it

A panel of a pairs plot IS a marginal, so
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)
exists exactly where the marginals do.

## Notation

\\A\\ is the retained subset of coordinates, \\A^c\\ its complement,
\\\mu\\ the location, \\\Sigma\\ the matrix the family carries, \\\nu\\
a Student t's degrees of freedom and \\\phi\\ a Dirichlet's
concentration.

## See also

[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
whose panels are these marginals, and
[`mv_marginal.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvGaussianDistrib.md),
[`mv_marginal.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvStudentTDistrib.md)
and
[`mv_marginal.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.multivariate_distrib.md)
for the methods and the refusal.

## Examples

``` r
d <- mvgaussian_distrib(3)
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
t3 <- mvstudent_t_distrib(3)
th <- as.list(stats::setNames(c(unlist(theta), 5), t3@params))
c(full = th$nu, marginal = mv_marginal(t3, th, c(1, 3))$theta$nu)
#>     full marginal 
#>        5        5 

# A Dirichlet's is a beta with the same concentration, and is a genuinely
# univariate object.
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
