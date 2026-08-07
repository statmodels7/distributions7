# No Skewness Without Saying Which One

Rejected. A scalar skewness for a vector response is not one quantity
but a choice among several – Mardia's, Malkovich-Afifi's, or the vector
of coordinatewise marginal skewnesses – and they do not agree. Returning
any of them under the bare name would be a wrong answer in the shape of
a right one, so the caller names the quantity it wants instead. Note
that
[`mv_marginal`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
is not a way round this for an elliptical family, whose marginal is a
smaller multivariate distribution and rejects in turn; it is for the
Dirichlet and the multinomial, whose marginals are univariate.

## Arguments

- x:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

Never returns; raises an error.

## See also

[`mv_marginal`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
