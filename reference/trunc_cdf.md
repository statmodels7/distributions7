# Distribution Function of a Truncated Distribution

\\F_T(q) = (F(q) - F(\ell^-))/Z\\, clamped to \\\[0, 1\]\\.

## Usage

``` r
trunc_cdf(distrib, q, theta, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- distrib:

  A truncated distribution object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log.p:

  Logical; whether to return the log probability.

## Value

A numeric vector.
