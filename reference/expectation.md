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
  [`distrib()`](https://statmodels7.github.io/distributions7/reference/distrib.md).

- f:

  The function whose expectation is taken, with signature
  `f(y, theta, ...)`. It is evaluated elementwise: `y` arrives as a
  numeric vector and every component of `theta`, like every argument
  passed through `...`, as a vector of the same length, so `f` must be
  vectorized in all of them jointly – which any expression built from
  arithmetic and the `distrib_*` generics already is.

- theta:

  A named list of parameters. Vectors are supported and are recycled
  against any vectors in `...`, so several parameter values can be
  handled in one call; all combinations share one batched evaluation.

- ...:

  Further arguments passed to `f`; their names must not clash with those
  of `theta`.

## Value

A numeric vector of expected values, one per parameter combination.

## See also

[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md),
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)

## Examples

``` r
expectation(poisson_distrib(), function(y, theta, k = 1) y^k,
            theta = list(mu = 10), k = 2)
#> [1] 110
```
