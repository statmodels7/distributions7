# The Pieces a Dirichlet Derivative Needs

Evaluates the mean vector, the concentration, the shapes \\\alpha =
\phi\mu\\ and the simplex's first two derivative arrays once, so that a
density or a derivative method computes them a single time and shares
them.

## Usage

``` r
dir_parts(distrib, theta)
```

## Arguments

- distrib:

  A
  [`DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/DirichletDistrib.md)
  object, read for its `param` and its parameter names.

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values, in the simplex's own order, followed by the concentration
  `phi`.

## Value

A named list with `mu`, the mean vector of length \\p\\ summing to one;
`phi`, the concentration; `alpha`, the shapes \\\phi\mu\\; `A`, the \\p
\times (p-1)\\ matrix of first derivatives of the mean in the free
values; `B`, the list of second-derivative vectors keyed by unordered
tuple; `idx`, those tuples; and `k`, the number of free mean values.

## Details

Two identities keep every formula on the family's pages short, and both
follow from differentiating \\\sum_j \mu_j = 1\\: once, the columns of
\\A = \partial\mu/\partial\eta\\ sum to zero; twice, so does every
second-derivative vector \\B\_{\cdot,kl}\\. They are what makes the
expected information closed form, since every term carrying the data is
a constant times one of those zero sums and drops out under expectation.

## See also

[`dir_b_index()`](https://statmodels7.github.io/distributions7/reference/dir_b_index.md)
for locating an entry of `B`,
[`distrib_gradient.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.DirichletDistrib.md)
for the first consumer, and
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
for the family.
