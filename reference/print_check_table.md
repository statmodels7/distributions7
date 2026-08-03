# Print a Validation Table

Renders what
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
found, in the same shape for a univariate and a multivariate
distribution.

## Usage

``` r
print_check_table(distrib, out, theta, n, nsim)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- out:

  The data frame of checks.

- theta:

  The parameters the checks were run at.

- n:

  The number of observations used.

- nsim:

  The Monte Carlo sample size used.

## Value

Invisibly `NULL`.

## See also

[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
