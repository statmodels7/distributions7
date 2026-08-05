# The Pieces a Dirichlet Derivative Needs

The mean vector, the concentration, the shapes \\\alpha = \phi\mu\\ and
the simplex's first two derivatives, computed once per call.

## Usage

``` r
dir_parts(distrib, theta)
```

## Arguments

- distrib:

  A
  [`DirichletDistrib`](https://statmodels7.github.io/distributions7/reference/DirichletDistrib.md)
  object.

- theta:

  A named list of parameters.

## Value

A list with `mu`, `phi`, `alpha`, `A` (a \\p \times (p-1)\\ matrix), `B`
(the second derivatives, keyed by tuple) and `idx` (their index tuples).

## Details

Two identities keep every formula short and are worth naming, both
following from \\\sum_j \mu_j = 1\\ differentiated once and twice: the
columns of \\A = \partial\mu/\partial\eta\\ sum to zero, and so does
every second-derivative vector. They are what make the expected
information closed form, since the terms carrying the data drop out
under expectation.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
