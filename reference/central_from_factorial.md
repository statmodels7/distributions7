# Central Moments From Falling Factorial Moments

Converts \\E\[Y^{(k)}\]\\ to the first four central moments, through the
raw moments \\m_1 = f_1\\, \\m_2 = f_2 + f_1\\, \\m_3 = f_3 + 3f_2 +
f_1\\ and \\m_4 = f_4 + 6f_3 + 7f_2 + f_1\\.

## Usage

``` r
central_from_factorial(f)
```

## Arguments

- f:

  A list of the four falling factorial moments.

## Value

A list with the mean and the second, third and fourth central moments.
