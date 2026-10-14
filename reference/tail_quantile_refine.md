# Refine a Quantile Far in a Tail

Corrects quantiles `x` of the requested tail where its log probability
is below `log_tail_cut`, by a safeguarded Newton iteration on \\r = \log
P(x) - \log p\\, with \\d\log P/dx = \pm f/P\\ from the family's own log
distribution function and log density.

## Usage

``` r
tail_quantile_refine(distrib, x, p, theta, lower.tail, log.p)
```

## Arguments

- distrib:

  A continuous distribution object.

- x:

  A numeric vector of starting quantiles.

- p, lower.tail, log.p:

  As in
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).

- theta:

  A named list of parameters.

## Value

A numeric vector of quantiles, `x` itself where the tail is not far.

## Details

On a positive support, and where \\x\\ lies on the side of the tail with
\\\|x\| \> 1\\, the variable is \\v = \log\|x\|\\, on which a polynomial
tail and a power of \\x\\ near zero are linear; elsewhere it is \\x\\.
Each step is the Newton step, limited to 10 units of \\v\\ (or to
\\\max(1, \|x\|)\\ in \\x\\), and replaced by the midpoint of the
current bracket of the root when it leaves it. A start that is not
finite or not inside the support is replaced by 1 on a positive support
and by 0 otherwise. The iteration stops when \\\|r\|\\ is within four
units in the last place of \\\log p\\; an element that does not reach it
in 200 steps is given to the default quantile of `continuous_distrib`.

## Examples

``` r
d <- invgauss1_distrib()
th <- list(mu = 1.7, phi = 2.3)
x0 <- suppressWarnings(statmod::qinvgauss(-800, mean = 1.7, dispersion = 2.3,
                                          lower.tail = FALSE, log.p = TRUE))
x <- distributions7:::tail_quantile_refine(d, x0, -800, th, FALSE, TRUE)
distrib_cdf(d, c(x0, x), th, lower.tail = FALSE, log.p = TRUE)
#> [1] -609.2349 -800.0000
```
