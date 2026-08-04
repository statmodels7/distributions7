# The Success Probability Behind a Geometric Mean

Converts the mean into the success probability the base R functions
take.

## Usage

``` r
geom_prob(mu)
```

## Arguments

- mu:

  The mean, a positive numeric vector.

## Value

A numeric vector of probabilities in \\(0, 1)\\.

## Details

The mean number of failures before the first success is \\(1-p)/p\\, so
\\p = 1/(1+\mu)\\. Writing it once keeps the four functions that call
stats from each repeating the algebra.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
