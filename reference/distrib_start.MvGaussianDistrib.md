# The Maximum Likelihood Estimate as a Starting Value

The sample mean and the sample covariance, which for an unstructured
covariance are the maximum likelihood estimate itself, so the fit begins
at the answer and confirms it in one step. For a structured covariance
they are the closest thing the matrix parameter can represent, which is
a good deal nearer than the origin.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  The response.

- n_start:

  Unused; one starting value is enough when it is the estimate.

- ...:

  Unused.

## Value

A list with one named parameter list.

## Details

When the matrix parameter parametrizes the precision the sample
covariance is inverted first, since that is the matrix the matrix
parameter has to represent.
