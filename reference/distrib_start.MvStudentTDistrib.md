# The Gaussian Estimate as a Starting Value for a t

The sample mean and the sample covariance, with the degrees of freedom
set where the family is heavy tailed but its second moment exists. The
gaussian estimate is the limit of this family as \\\nu\\ grows, so it is
the right place to start looking for a finite one.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  The response.

- n_start:

  Unused.

- ...:

  Unused.

## Value

A list with one named parameter list.

## Details

The scale matrix is the covariance divided by \\\nu/(\nu-2)\\, and that
factor is applied, since a starting value that confused the two would
begin with a scale a third too large.
