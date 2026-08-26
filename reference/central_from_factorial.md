# Central Moments From Falling Factorial Moments

Converts the first four falling factorial moments into the mean and the
second, third and fourth central moments. The raw moments come first,
through the Stirling numbers of the second kind written out: \$\$m_1 =
f_1,\quad m_2 = f_2 + f_1,\quad m_3 = f_3 + 3f_2 + f_1,\quad m_4 = f_4 +
6f_3 + 7f_2 + f_1,\$\$ and the central moments follow by the usual
binomial expansion about \\m_1\\.

## Usage

``` r
central_from_factorial(f)
```

## Arguments

- f:

  An unnamed list of four, the falling factorial moments of order 1 to 4
  in that order, as
  [`betabinom_factorial_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_factorial_moments.md)
  returns. They are read by position.

## Value

A named list with `mean` and the central moments `c2`, `c3` and `c4`,
each a numeric vector the length of the inputs.

## See also

[`betabinom_factorial_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_factorial_moments.md)
for the input;
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md),
which chains the two.

## Examples

``` r
f <- distributions7:::betabinom_factorial_moments(2, 2, 10)
distributions7:::central_from_factorial(f)$mean
#> [1] 5
```
