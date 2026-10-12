# Tail Probabilities on the Log Scale

Three internal helpers for the distribution and quantile functions of
the families with a closed form. `log1mexp_neg()` computes \\\log(1 -
e^a)\\ for \\a \le 0\\, through `log(-expm1(a))` when \\a \> -\log 2\\
and through `log1p(-exp(a))` otherwise, so that the result keeps its
relative accuracy at both ends (Maechler, 2012). `cdf_from_log_tails()`
returns the probability that `lower.tail` and `log.p` select from the
logarithms of the two tails, both computed by the caller.
`quantile_log_tails()` converts a probability given with `lower.tail`
and `log.p` into the logarithms of the lower and the upper tail, without
forming \\1 - p\\ or exponentiating a logarithm; a probability outside
\\\[0, 1\]\\ gives `NaN` for both.

## Usage

``` r
log1mexp_neg(a)

cdf_from_log_tails(log_lower, log_upper, lower.tail, log.p)

quantile_log_tails(p, lower.tail, log.p)
```

## Arguments

- a:

  A numeric vector of non-positive numbers.

- log_lower, log_upper:

  Numeric vectors, the logarithms of \\P(Y \le q)\\ and \\P(Y \> q)\\.

- lower.tail, log.p:

  Logical of length 1, with the meaning they have in
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
  and
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).

- p:

  A numeric vector of probabilities, or of their logarithms when
  `log.p = TRUE`.

## Value

`log1mexp_neg()` and `cdf_from_log_tails()` return a numeric vector of
the length of their input; `quantile_log_tails()` returns a list with
the numeric vectors `lower` and `upper`.

## References

Maechler, M. (2012). Accurately computing \\\log(1 - \exp(-\|a\|))\\.
Vignette of the R package Rmpfr.

## Examples

``` r
# The helper keeps the digits that the direct formula loses (it gives -Inf).
distributions7:::log1mexp_neg(-1e-20)
#> [1] -46.0517
log(1 - exp(-1e-20))
#> [1] -Inf

# The upper tail of a probability given on the log scale.
distributions7:::quantile_log_tails(log(1e-30), lower.tail = FALSE,
                                    log.p = TRUE)
#> $lower
#> [1] -1e-30
#> 
#> $upper
#> [1] -69.07755
#> 
```
