# The Pieces a Multinomial Derivative Needs

Evaluates the probability vector and the simplex's first two derivative
arrays once, so that a mass function or a derivative method computes
them a single time and shares them.

## Usage

``` r
mn_parts(distrib, theta)
```

## Arguments

- distrib:

  A
  [`MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/MultinomialDistrib.md)
  object, read for its `param`.

- theta:

  A named list of parameters on the parameter scale, the simplex's free
  values in its own order.

## Value

A named list with `prob`, the probability vector of length \\p\\ summing
to one; `A`, the \\p \times (p-1)\\ matrix of first derivatives in the
free values; `B`, the list of second-derivative vectors keyed by
unordered tuple; `idx`, those tuples; and `k`, the number of free
values.

## Details

The probabilities sum to one at every free vector, so differentiating
that identity shows that the columns of \\A = \partial p/\partial\eta\\
sum to zero and so does every second-derivative vector
\\B\_{\cdot,kl}\\. The second of those makes the expected information
closed form: under expectation the term carrying \\B\\ becomes \\n\sum_j
B\_{j,kl}\\, which is zero.

## See also

[`distrib_gradient.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MultinomialDistrib.md)
for the first consumer,
[`dir_parts()`](https://statmodels7.github.io/distributions7/reference/dir_parts.md)
for the Dirichlet's counterpart, and
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
for the family.
