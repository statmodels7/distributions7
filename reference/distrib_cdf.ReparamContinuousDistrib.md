# Distribution Function of a Reparametrized Distribution

The parent's distribution function at the mapped parameters.

## Usage

``` r
reparam_cdf(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- distrib:

  A reparametrized distribution.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters.

- lower.tail:

  Logical; if `TRUE`, probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, returns log-probabilities.

- ...:

  Passed to the parent.

## Value

A numeric vector.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
