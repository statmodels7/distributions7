# Quantile Function of a Reparametrized Distribution

The parent's quantile function at the mapped parameters.

## Usage

``` r
reparam_quantile(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...)
```

## Arguments

- distrib:

  A reparametrized distribution.

- p:

  A numeric vector of probabilities.

- theta:

  A named list of the new parameters.

- lower.tail:

  Logical; if `TRUE`, probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is a log-probability.

- ...:

  Passed to the parent.

## Value

A numeric vector.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
