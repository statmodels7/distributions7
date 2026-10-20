# The Support Points of a Discrete Multivariate Distribution

Returns the points on which a discrete multivariate distribution places
mass, as a matrix with one row per point. With the support enumerated,
every expectation is an exact sum, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
verifies the total mass by addition instead of by sampling.

## Usage

``` r
mv_support(distrib, theta, ...)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list or vector of parameters. The shipped method does not read
  it, because the support of the multinomial does not depend on them.

- ...:

  Passed to methods. The shipped method does not read it.

## Value

A numeric matrix with one row per support point and one column per
coordinate.

## Why the generic exists

A univariate discrete distribution does not need it, since its support
is a range of integers, which the package traverses. The support of a
vector is a set whose shape depends on the family, so each discrete
multivariate family enumerates its own.

## The multinomial's support

For \\p\\ coordinates and \\n\\ trials it is the set of weak
compositions of \\n\\ into \\p\\ parts, \$\$\mathcal{S} = \Bigl\\y \in
\mathbb{N}\_0^{p} : \textstyle\sum\_{j=1}^{p} y_j = n\Bigr\\, \qquad
\lvert\mathcal{S}\rvert = \binom{n + p - 1}{p - 1},\$\$ enumerated by
[`numericals7::compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.html).
Every expectation is then the finite sum \\\sum\_{y \in \mathcal{S}}
g(y) f(y; \theta)\\. The count grows quickly: 15 points at \\n = 4, p =
3\\, and 5151 at \\n = 100, p = 3\\.

## The method of the class

It signals an error. A continuous family has no such set, and a discrete
family whose support is infinite has no finite matrix to return.
[`has_mv_support()`](https://statmodels7.github.io/distributions7/reference/has_mv_support.md)
tests whether a family has a method of its own.

## Notation

\\\mathcal{S}\\ is the support, \\p\\ the dimension, \\n\\ the number of
trials of a multinomial and \\f\\ the mass function.

## See also

[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md),
the one family with a method,
[`has_mv_support()`](https://statmodels7.github.io/distributions7/reference/has_mv_support.md)
for the test,
[`numericals7::compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.html)
for the enumeration, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
for the consumer.

## Examples

``` r
d <- multinomial_distrib(3, size = 4)
theta <- list(probs_alr1 = 0.4, probs_alr2 = -0.3)
sup <- mv_support(d, theta)

# Each row is one point, and each row sums to the number of trials.
dim(sup)
#> [1] 15  3
unique(rowSums(sup))
#> [1] 4
head(sup, 4)
#>      [,1] [,2] [,3]
#> [1,]    0    0    4
#> [2,]    0    1    3
#> [3,]    0    2    2
#> [4,]    0    3    1

# The count is the number of weak compositions.
c(got = nrow(sup), expected = choose(4 + 3 - 1, 3 - 1))
#>      got expected 
#>       15       15 

# The mass over the support sums to one exactly.
sum(distrib_pdf(d, sup, theta))
#> [1] 1

# A continuous family has no such set, and the call signals an error.
try(mv_support(mvgaussian1_distrib(2),
               list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
                    sigma_log_L2 = 0, sigma_L2.1 = 0)))
#> Error : 'multivariate gaussian [2d, sigma=log_cholesky]' does not enumerate a support. A continuous family has no such set,
#>   and a discrete one whose support is infinite has no finite matrix to
#>   return; a family that does have one registers this generic.
```
