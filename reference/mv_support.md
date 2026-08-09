# The Support Points of a Discrete Multivariate Distribution

The points a discrete multivariate distribution places mass on, as a
matrix with one row per point.

## Usage

``` r
mv_support(distrib, theta, ...)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters. Families whose support does not depend on
  them ignore it.

- ...:

  Passed to methods.

## Value

A matrix with one row per support point and one column per coordinate.

## Details

A univariate discrete distribution needs no such generic: its support is
a stretch of the integers and the package walks it. On a vector the
support is a set whose shape depends on the family — the multinomial's
is the compositions of its size — and enumerating it is what lets an
expectation be an exact sum and the validator check the total mass by
addition rather than by sampling.

The base class rejects. A continuous family has no such set, and a
discrete one whose support is infinite has no finite matrix to return;
either way an answer would be a fiction, and the caller is better told.

For the multinomial on \\p\\ coordinates with \\n\\ trials the support
is the weak compositions of \\n\\ into \\p\\ parts,

\$\$\mathcal{S} = \Bigl\\y \in \mathbb{N}\_0^{p} :
\textstyle\sum\_{j=1}^{p} y_j = n\Bigr\\, \qquad \lvert\mathcal{S}\rvert
= \binom{n + p - 1}{p - 1},\$\$

enumerated by
[`compositions`](https://statmodels7.github.io/numericals7/reference/compositions.html).
Every expectation is then the finite sum \\\sum\_{y \in \mathcal{S}}
g(y) f(y; \theta)\\.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md),
[`compositions`](https://statmodels7.github.io/numericals7/reference/compositions.html)

## Examples

``` r
d <- multinomial_distrib(3, size = 4)
nrow(mv_support(d, list(probs_alr1 = 0, probs_alr2 = 0)))
#> [1] 15
```
