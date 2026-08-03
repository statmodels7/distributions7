# Expected Value of a Function of a Random Variable

Computes \\E\[f(Y)\]\\ under the distribution: by numerical integration
for a continuous distribution and by series summation for a discrete
one, with the methods described on their own pages.

## Usage

``` r
expectation(distrib, f, theta, ...)
```

## Arguments

- distrib:

  An object inheriting from
  [`distrib`](https://statmodels7.github.io/distributions7/reference/distrib.md).

- f:

  The function whose expectation is taken, with signature
  `f(y, theta, ...)`.

- theta:

  A named list of parameters. Vectors are supported and are recycled
  against any vectors in `...`, so several parameter values can be
  handled in one call.

- ...:

  Further arguments passed to `f`; their names must not clash with those
  of `theta`.

## Value

A numeric vector of expected values, one per parameter combination.

## Examples

``` r
expectation(poisson_distrib(), function(y, theta, k = 1) y^k,
            theta = list(mu = 10), k = 2)
#> [1] 110
```
