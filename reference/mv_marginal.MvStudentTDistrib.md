# Marginal of a Multivariate Student t

Returns the marginal law of a subset of coordinates, which for this
family is again a Student t: the subvector of the location, the
corresponding block of the scale matrix, and the same degrees of
freedom. The mixing variable of the scale-mixture representation is
shared by every coordinate, so it survives integrating any of them out.
The returned object is a new
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
of the reduced dimension when `distrib` is a
[MvStudentT1Distrib](https://statmodels7.github.io/distributions7/reference/MvStudentT1Distrib.md),
and a new
[`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
when it is a
[MvStudentT2Distrib](https://statmodels7.github.io/distributions7/reference/MvStudentT1Distrib.md),
both on an unstructured log-Cholesky matrix. The free values are
computed by
[`parameters7::param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.html)
from the block of the scale matrix (from its inverse for the second
family), so they are not a subset of the original's.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
  or
  [`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- theta:

  A named list of parameters, each component a single number.

- which:

  An integer vector of distinct coordinates to keep, each between 1 and
  \\p\\, already checked by the generic
  [`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md).
  The marginal follows the order given.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `distrib`, an `MvStudentTDistrib` of dimension
`length(which)` of the parent's subclass, and `theta`, its parameters as
a named list.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom and \\p\\ the dimension of the full law.

## See also

[`mv_sigma.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
for the matrix it takes a block of,
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
whose diagonal panels are these marginals, and
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(3)
theta <- as.list(stats::setNames(
  c(1, -2, 0.5, 0.1, -0.2, 0.3, 0.4, -0.1, 0.2, 6), d@params))

m <- mv_marginal(d, theta, c(1, 3))
m$distrib@n_dim
#> [1] 2

# The degrees of freedom are unchanged, and the scale block is the
# submatrix of the full one.
c(full = theta$nu, marginal = m$theta$nu)
#>     full marginal 
#>        6        6 
all.equal(mv_sigma(m$distrib, m$theta), mv_sigma(d, theta)[c(1, 3), c(1, 3)],
          check.attributes = FALSE)
#> [1] TRUE

# A single coordinate is the univariate t, against stats::dt scaled.
m1 <- mv_marginal(d, theta, 1)
s1 <- sqrt(mv_sigma(d, theta)[1, 1])
c(ours = distrib_pdf(m1$distrib, 2.4, m1$theta, log = TRUE),
  dt = dt((2.4 - 1) / s1, df = 6, log = TRUE) - log(s1))
#>      ours        dt 
#> -1.889948 -1.889948 
```
