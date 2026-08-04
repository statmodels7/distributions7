# Scale Standard Deviations and Correlations of a Multivariate t

The square roots of the diagonal of the **scale** matrix and the
correlations it implies. The correlations are those of the response as
well, since the covariance is \\\nu\Sigma/(\nu-2)\\ and a positive
multiple of a matrix leaves its correlations alone; the diagonal
quantities are not standard deviations of the response and are named to
say so.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

A list as described in
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).
