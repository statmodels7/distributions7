# Falling Factorial Moments of a Beta-Binomial

Returns the first four falling factorial moments \$\$E\[Y^{(k)}\] =
n^{(k)} \prod\_{j=0}^{k-1} \frac{a+j}{a+b+j}, \qquad k = 1, \ldots,
4,\$\$ with \\x^{(k)} = x(x-1)\cdots(x-k+1)\\. These are the quantities
a beta-binomial has in closed form; the raw and central moments follow
from them through
[`central_from_factorial()`](https://statmodels7.github.io/distributions7/reference/central_from_factorial.md).

## Usage

``` r
betabinom_factorial_moments(a, b, n)
```

## Arguments

- a:

  The first shape of the mixing beta, a single positive number.

- b:

  The second shape, a single positive number.

- n:

  The number of trials, a single non-negative whole number. Factorial
  moments of order above `n` are zero, the falling factorial \\n^{(k)}\\
  vanishing there.

## Value

An **unnamed** list of four numbers, the falling factorial moments of
order 1 to 4 in that order, reached by position.

## Notation

\\a \> 0\\ and \\b \> 0\\ are the two shapes of the mixing beta, \\n\\
the number of trials, and \\x^{(k)}\\ the falling factorial.

## See also

[`central_from_factorial()`](https://statmodels7.github.io/distributions7/reference/central_from_factorial.md),
which converts these;
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md),
the wrapper the moment methods call.

## Examples

``` r
# At a = b the mixing beta is symmetric, so the mean is half the trials.
distributions7:::betabinom_factorial_moments(2, 2, 10)[[1]]
#> [1] 5
```
