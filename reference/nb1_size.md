# The Size Behind an NB1 Mean

The number of successes \\r = \mu/\theta\\ that the base R negative
binomial functions take.

## Usage

``` r
nb1_size(mu, theta)
```

## Arguments

- mu:

  The mean, a positive numeric vector.

- theta:

  The dispersion, a positive numeric vector.

## Value

A numeric vector.

## Details

Requiring the variance to be \\\mu(1+\theta)\\ fixes the success
probability at \\1/(1+\theta)\\, and the mean then determines the size.
It is this that puts \\\mu\\ inside the gamma functions and makes the
family different from the quadratic form rather than a reparametrization
of it.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md),
[`nb1_prob`](https://statmodels7.github.io/distributions7/reference/nb1_prob.md)
