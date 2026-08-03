# Quantile Function of a Truncated Distribution

Inverts
[`trunc_cdf`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
through the parent's quantile function.

## Usage

``` r
trunc_quantile(distrib, p, theta, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- distrib:

  A truncated distribution object.

- p:

  A numeric vector of probabilities.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; whether `p` is a lower-tail probability.

- log.p:

  Logical; whether `p` is given on the log scale.

## Value

A numeric vector of quantiles.

## Details

Inverse transform on the parent: \\F_T(q) = p\\ exactly when \\F(q) =
F(\ell^-) + pZ\\, so no root-finding of its own is needed. The
generalized inverse of a discrete cdf satisfies the same relation, so
the discrete case needs no separate treatment.
